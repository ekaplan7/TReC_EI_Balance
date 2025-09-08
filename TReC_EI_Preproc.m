%code adapted from existing preproc pipeline from Dr. Danielle Sliva 
clear all;
%eeglab

%eeglab specific format 
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

rawdatadir = '/Users/elizabethkaplan/Desktop/Greene/source_data';
subject_dir = '/Users/elizabethkaplan/Desktop/Greene/subject_folders/TReC_12237';
analysis_dir = '/Users/elizabethkaplan/Desktop/Greene/analysis_files';

%% Load in data - .set file 

setfile = 'TReC_12237_PreEEG.set';   

% load the dataset into EEGLAB format 
EEG = pop_loadset('filename', setfile, 'filepath', subject_dir);

%save in this format for later visual inspection
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 1);

%% Load in data - .cnt file

cnt_file = fullfile(subject_dir, 'study-353_sub-12237_task-eeg_ses-post.cnt');
EEG = pop_fileio(cnt_file);

%% add channel and block info -

chanloc_file = fullfile(fileparts(which('eeglab.m')), 'sample_locs', 'standard_waveguard64.elc');
EEG = pop_chanedit(EEG, 'load', {chanloc_file, 'filetype', 'autodetect'});
EEG_mask = EEG;

%update block times for each subj!! 
Blocks_onset_ms   = [60298 250298 438298 626298 814298 1002298];
Blocks_offset_ms  = [240298 430298 618298 806298 994298 1182298];

%% 
function EEG = insert_block_events(EEG, blocks_onset_ms, blocks_offset_ms)
% insert_block_events Add trial onset/offset events to EEG struct
% Inputs:
%   EEG              - EEGLAB dataset struct, already loaded
%   blocks_onset_ms  - vector of onset times (ms)
%   blocks_offset_ms - vector of offset times (ms)
% vectors must be same length
%
% Output:
%   EEG              - EEGLAB dataset struct with added trial events

fs = EEG.srate;
video_delay_sec = 0.8960; 
delay_samples = video_delay_sec * fs; 

onset_samples = (blocks_onset_ms / 1000) * fs + delay_samples;
offset_samples = (blocks_offset_ms / 1000) * fs + delay_samples;

% insert onsets
for idx = 1:length(onset_samples)
    latency_samps = onset_samples(idx);
    EEG.event(end+1).type = 'block_onset';
    EEG.event(end).latency = latency_samps;
    EEG.event(end).urevent = length(EEG.urevent) + 1;
end

% insert offsets
for idx = 1:length(offset_samples)
    latency_samps = offset_samples(idx);
    EEG.event(end+1).type = 'block_offset';
    EEG.event(end).latency = latency_samps;
    EEG.event(end).urevent = length(EEG.urevent) + 1;
end

EEG = eeg_checkset(EEG, 'eventconsistency');

end

EEG = insert_block_events(EEG, Blocks_onset_ms, Blocks_offset_ms);

%% Make sure trial events were correctly loaded to EEG structure 

EEG.event.type

for i = 1:min(10, length(EEG.event))
    fprintf('Event %d: type=%s, latency=%.1f\n', ...
        i, EEG.event(i).type, EEG.event(i).latency);
end

%% Load tic file - containing ALL tics 

tic_event_file = '/Users/elizabethkaplan/Desktop/Greene/subject_folders/TReC_12237/TReC_12237_PreEEG_all_timestamps_sorted2.txt';

if ~isfile(tic_event_file)
    error('Timestamp file not found for subject %s, timepoint %s', subjectID, timepoint);
end

tic_event_table = readtable(tic_event_file);
tic_event_table.type = string(tic_event_table.type); 

% Convert to seconds 
tic_event_table.latency_sec = tic_event_table.latency * 60;

% tic types
motor_onset_codes = {'101','102','103','104','105','106','107','108','109','110','111','112'};
motor_offset_codes = {'201','202','203','204','205','206','207','208','209','210','211','113'};
vocal_onset_codes = {'301','302','303','304','305'};
vocal_offset_codes = {'401','402','403','404','405'};

tic_periods = [];

function periods = pair_onset_offset(onsets, offsets)
    periods = [];
    while ~isempty(onsets)
        onset_time = onsets.latency_sec(1);
        % Find first offset after this onset
        offset_idx = find(offsets.latency_sec > onset_time, 1, 'first');
        if isempty(offset_idx)
            warning('No matching offset found for onset at %.2f sec', onset_time);
            onsets(1,:) = [];
            continue;
        end
        offset_time = offsets.latency_sec(offset_idx);
        periods = [periods; onset_time, offset_time];
        % Remove used rows
        onsets(1,:) = [];
        offsets(offset_idx,:) = [];
    end
end

%% Remove Tics 
motor_onsets = tic_event_table(ismember(tic_event_table.type, motor_onset_codes), :);
motor_offsets = tic_event_table(ismember(tic_event_table.type, motor_offset_codes), :);
vocal_onsets = tic_event_table(ismember(tic_event_table.type, vocal_onset_codes), :);
vocal_offsets = tic_event_table(ismember(tic_event_table.type, vocal_offset_codes), :);

tic_periods = [tic_periods; pair_onset_offset(motor_onsets, motor_offsets)];
tic_periods = [tic_periods; pair_onset_offset(vocal_onsets, vocal_offsets)];

% Merge overlapping periods
if ~isempty(tic_periods)
    tic_periods = sortrows(tic_periods, 1);
    merged_periods = tic_periods(1,:);
    for i = 2:size(tic_periods,1)
        if tic_periods(i,1) <= merged_periods(end,2)  % overlap
            merged_periods(end,2) = max(merged_periods(end,2), tic_periods(i,2));
        else
            merged_periods = [merged_periods; tic_periods(i,:)];
        end
    end
    tic_periods = merged_periods;
end

%actually remove 
if ~isempty(tic_periods)
    EEG = pop_select(EEG, 'notime', tic_periods);  % EEGLAB expects Nx2 transposed
    EEG = eeg_checkset(EEG);
    fprintf('Removed %d tic periods from EEG.\n', size(tic_periods,1));
else
    fprintf('No tic periods found to remove.\n');
end

%% 3. Remove DC offset (channel mean) 

for iChan = 1:size(EEG.data, 1)
    EEG.data(iChan, :) = single(double(EEG.data(iChan, :)) - ...
                                mean(double(EEG.data(iChan, :))));
end

%% 4. Notch filter 60 Hz - to remove line noise (noise from electricity at 60 Hz)

EEG_params.notch.low_cutoff = 58; 
EEG_params.notch.high_cutoff = 62;
notch_filter = 1;

EEG = pop_eegfiltnew(EEG, EEG_params.notch.low_cutoff, ...
    EEG_params.notch.high_cutoff, [], notch_filter);
EEG = eeg_checkset(EEG);

% Save updated EEG_params
%save(params_file, 'EEG_params');

%% 5. High and low pass filtering - adjust as needed

EEG = pop_basicfilter(EEG, 1:EEG.nbchan, ...
        'Cutoff' , .1, ...
        'Design' , 'butter' , ...
        'Filter' , 'highpass', ...
        'Order'  , 2 , ...
        'RemoveDC','on', ...
        'Boundary','none');

EEG = pop_basicfilter(EEG, 1:EEG.nbchan, ...
        'Cutoff' , 70 , ...
        'Design'  , 'butter' , ...
        'Filter'  , 'lowpass', ...
        'Order'   , 4 , ...
        'Boundary', 'none');

%% Run sliding window  

fs = EEG.srate;                    
winSize_sec = 0.05;                
step_sec = 0.025;                  
threshold = 200;                    
channels = 1:EEG.nbchan;          

winSamples = round(winSize_sec * fs);
stepSamples = round(step_sec * fs);
nSamples = EEG.pnts;

rejMask = false(1, nSamples);

% Sliding window
for startS = 1:stepSamples:(nSamples-winSamples+1)
    idx = startS:(startS+winSamples-1);
    if any(abs(EEG.data(channels, idx)) > threshold, 'all')
        rejMask(idx) = true;  % Mark samples for rejection
    end
end

% convert to secs
diffMask = diff([0, rejMask, 0]);
starts = find(diffMask == 1);
ends   = find(diffMask == -1) - 1;
rejectedRanges_sec = [starts; ends]' / fs;
rejectedRanges_samples = round(rejectedRanges_sec * fs);

% Remove
EEG_clean = EEG;
EEG_clean = pop_select(EEG_clean, 'notime', rejectedRanges_sec);

EEG_clean = eeg_checkset(EEG_clean);

%% Remove bad channels 
%specifically looking for channels that do not follow pattern of power
%spectra pl
figure; pop_spectopo(EEG_clean, 1, [], 'EEG' , 'percent',15,'freq',...
        [10 20 30],'freqrange',[1 70],'electrodes','off');

% Signal
%pop_eegplot( EEG_clean, 1, 1, 1);

disp('----------------------')
disp('----------------------')
disp('Manual Input Required')
disp('----------------------')
disp('----------------------')
bad = input('Select bad channels:');

EEG=pop_select(EEG_clean, 'nochannel',bad);

%% ICA to remove blinks and muscle artifact 

%create copy of the data highpassed at 1 Hz - improves identification of
%eye blinks and extreme muscle artifacts (aka smooths out signal so these are
%more obvious)
EEG = pop_basicfilter(EEG,1:EEG.nbchan,'Cutoff',1,'Filter','highpass','Design','butter','Order',2);

% Down‑sample to 256 Hz first (faster)
EEG = pop_resample(EEG,256);

%last input is the num of ICs you want computed -- this number must be less
%than the number of channels you actually have so may need updating
%depending on how many channels removed 
EEG = pop_runica(EEG,'icatype','picard','pca',60);

%longer ICA version - takes alot more time to run
%EEG_ica = pop_runica(EEG_ica,'extended',1,'icatype','runica','pca',EEG_ica.nbchan);

%convert weights back to original dataframe
%EEG.icaweights = EEG.icaweights;
%EEG.icasphere  = EEG_ica.icasphere;
%EEG.icachansind  = EEG_ica.icachansind;
%EEG.icawinv = EEG_ica.icawinv;
%EEG = eeg_checkset(EEG);

%% Visualize ICs 

figure;
pop_selectcomps(EEG, 1:35);  % Show the first 35 components (you can adjust this range if needed, but components larger than this are more unreliable)
pop_eegplot(EEG, 0, 1, 1);    

disp('----------------------')
disp('----------------------')
disp('Manual Input Required: Select Components to Reject')
disp('----------------------')
disp('----------------------')

% Prompt for user input to select components to reject
% You can select multiple components using a vector (e.g., [1 2 5 8])
reject_components = input('Enter the list of components to reject (e.g., [1 2 5 8]): ');

% Print the list of rejected components to the output for reference
%txt = sprintf('Reject Components: %s\n', num2str(reject_components));
%fprintf(out_txt, '%s\n ', txt); 
%fprintf('%s', txt);

% Remove the selected components
EEG = pop_subcomp(EEG, reject_components, 0); 

%% 13. Interpolate the rejected channels using spherical interpolation 

EEG = pop_interp(EEG, bad, 'spherical'); 

%% 14. Visual Inspection of data - "final pass" 
%can manually remove bad sections of data - be sure to hit "reject" in the bottom right corner for EEGLAB to
%automatically remove marked segments

EEG = eeg_checkset(EEG);

pop_eegplot(EEG, 1, 1, 1);
eegh

%% Save processed data with chanloc info to local device

EEG = pop_saveset(EEG, 'filename', 'TRec_12237_EI', 'filepath', '/Users/elizabethkaplan/Desktop/Greene/subject_folders/TReC_12237');
 