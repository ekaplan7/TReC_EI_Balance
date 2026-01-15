%code adapted from existing preproc pipeline from Dr. Danielle Sliva 
clear all;
%eeglab

%eeglab specific format 
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%UPDATE FILE PATH 
rawdatadir = '/Users/elizabethkaplan/Desktop/Greene/source_data';
subject_dir = '/Users/elizabethkaplan/Desktop/Greene/subject_folders/TReC_12297';
analysis_dir = '/Users/elizabethkaplan/Desktop/Greene/analysis_files';

%% Load in data - .set file 

%UPDATE FILE PATH 
setfile = 'TReC_12297_preEEG.set';   

% load the dataset into EEGLAB format 
EEG = pop_loadset('filename', setfile, 'filepath', subject_dir);

%save in this format for later visual inspection
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 1);

%% Load in data - .cnt file

%cnt_file = fullfile(subject_dir, 'study-353_sub-12237_task-eeg_ses-post.cnt');
%EEG = pop_fileio(cnt_file);

%% add channel and block info -

chanloc_file = fullfile(fileparts(which('eeglab.m')), 'sample_locs', 'standard_waveguard64.elc');
EEG = pop_chanedit(EEG, 'load', {chanloc_file, 'filetype', 'autodetect'});
EEG_mask = EEG;

%update block times for each subj!! -- in the primarycoded.csv file from
%Box
%Blocks_onset_min   = [1.1255 4.323 7.4297 10.5857 13.7045 16.8413];

%Blocks_offset_min  = [4.2183 7.4158 10.5244 13.6730 16.7947 19.9329];

%Blocks_onset_ms  = Blocks_onset_min  * 60000;
%Blocks_offset_ms = Blocks_offset_min * 60000;

%UPDATE NUMBERS 
Blocks_onset_ms = [106000 296000 485000 673000 861000 1049000];
Blocks_offset_ms = [286100 476000 665000 853000 1041000 1229000];

%calculate actual timing, accounting for offsets
%UPDATE NUMBERS
fs = EEG.srate;
video_delay_sec = 6.528; 
delay_samples = video_delay_sec * fs; 

% Convert ms → seconds
block_onsets_sec  = Blocks_onset_ms  / 1000;
block_offsets_sec = Blocks_offset_ms / 1000;

% Add video delay (still in seconds)
block_onsets_sec  = block_onsets_sec  + video_delay_sec;
block_offsets_sec = block_offsets_sec + video_delay_sec;

% Convert seconds → samples
onset_samples  = block_onsets_sec  * fs;
offset_samples = block_offsets_sec * fs;

% Make sure they are row vectors
onset_samples  = onset_samples(:)';
offset_samples = offset_samples(:)';

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

%% 5. High and low pass filtering

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

%% Load tic file - containing ALL tics 

%UPDATE FILE PATH 
tic_event_file = '/Users/elizabethkaplan/Desktop/Greene/subject_folders/TReC_12454/TReC_12454_PostEEG_all_timestamps_sorted.txt';

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

%% Load in tics 
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

tic_periods = tic_periods;

%% Run sliding window  

fs = EEG.srate;
winSize_sec = 0.05;
step_sec = 0.025;
threshold = 200;
channels = 1:EEG.nbchan;

winSamples  = round(winSize_sec * fs);
stepSamples = round(step_sec * fs);
nSamples    = EEG.pnts;

rejMask = false(1, nSamples);

for startS = 1:stepSamples:(nSamples - winSamples + 1)
    idx = startS:(startS + winSamples - 1);
    if any(abs(EEG.data(channels, idx)) > threshold, 'all')
        rejMask(idx) = true;
    end
end

% padding
pad_sec = 0.05; 
padSamp = round(pad_sec * fs);
if padSamp > 0
    rejMask = conv(double(rejMask), ones(1, 2*padSamp+1), 'same') > 0;
end

% convert to sample ranges 
diffMask = diff([0, rejMask, 0]);
starts = find(diffMask == 1);
ends   = find(diffMask == -1) - 1;

% Convert to seconds
rejectedRanges_sec = [starts(:) ends(:)] / fs;

%% Combine tic ranges and artifact ranges 

all_remove = [tic_periods; rejectedRanges_sec];
all_remove = sortrows(all_remove, 1);

% Merge overlaps so they don’t double-count
merged = all_remove(1,:);
for i = 2:size(all_remove,1)
    if all_remove(i,1) <= merged(end,2)
        merged(end,2) = max(merged(end,2), all_remove(i,2));
    else
        merged = [merged; all_remove(i,:)];
    end
end
all_remove = merged;

block_onsets_sec  = onset_samples  / fs;
block_offsets_sec = offset_samples / fs;

%check whether removal ranges overlap with block timestamps 
fprintf('\n=== Checking block onset/offset markers against removed segments ===\n\n');

protect_mask_onset  = true(size(block_onsets_sec));
protect_mask_offset = true(size(block_offsets_sec));

for j = 1:size(all_remove,1)
    rem_start = all_remove(j,1);
    rem_end   = all_remove(j,2);

    % Flag any onset that falls inside this removal window
    inside_onset  = block_onsets_sec >= rem_start & block_onsets_sec <= rem_end;
    protect_mask_onset(inside_onset) = false;

    % Flag any offset that falls inside this removal window
    inside_offset = block_offsets_sec >= rem_start & block_offsets_sec <= rem_end;
    protect_mask_offset(inside_offset) = false;
end

for i = 1:length(block_onsets_sec)
    if ~protect_mask_onset(i)
        fprintf('⚠️  Block %d onset (%.2f s) falls INSIDE a removed segment!\n', ...
            i, block_onsets_sec(i));
    else
        fprintf('✅ Block %d onset (%.2f s) is safe.\n', i, block_onsets_sec(i));
    end
end

for i = 1:length(block_offsets_sec)
    if ~protect_mask_offset(i)
        fprintf('⚠️  Block %d offset (%.2f s) falls INSIDE a removed segment!\n', ...
            i, block_offsets_sec(i));
    else
        fprintf('✅ Block %d offset (%.2f s) is safe.\n', i, block_offsets_sec(i));
    end
end

%% protect onset and offset markers 

protect_margin = 0.004; % 4 ms margin to ensure no overlap - may need to tweak based on subjs tic behavior 

% Make a copy so comparisons use the original windows
original_remove = all_remove;

for i = 1:length(block_onsets_sec)
    for j = 1:size(original_remove, 1)

        % --- Protect onset ---
        if block_onsets_sec(i) >= original_remove(j,1) && block_onsets_sec(i) <= original_remove(j,2)
            fprintf('⚠️ Onset %.3f s inside removed segment (%.3f–%.3f s). Trimming.\n', ...
                block_onsets_sec(i), original_remove(j,1), original_remove(j,2));

            % Shorten removal window so it ends just before the onset
            all_remove(j,2) = block_onsets_sec(i) - protect_margin;
        end

        % --- Protect offset ---
        if block_offsets_sec(i) >= original_remove(j,1) && block_offsets_sec(i) <= original_remove(j,2)
            fprintf('⚠️ Offset %.3f s inside removed segment (%.3f–%.3f s). Trimming.\n', ...
                block_offsets_sec(i), original_remove(j,1), original_remove(j,2));

            % Shorten removal window so it starts just after the offset
            all_remove(j,1) = block_offsets_sec(i) + protect_margin;
        end

    end
end

%% 
function EEG = insert_block_events(EEG, onset_samples, offset_samples)
% insert_block_events Add trial onset/offset events to EEG struct
% Inputs:
%   EEG             - EEGLAB dataset struct, already loaded
%   onset_samples   - vector of block onset latencies (in samples)
%   offset_samples  - vector of block offset latencies (in samples)
%
% Output:
%   EEG             - EEGLAB dataset struct with added events

if length(onset_samples) ~= length(offset_samples)
    error('Onset and offset vectors must be the same length.');
end

% --- Insert onsets ---
for idx = 1:length(onset_samples)
    EEG.event(end+1).type     = 'block_onset';
    EEG.event(end).latency    = onset_samples(idx);
    EEG.event(end).urevent    = length(EEG.event) + 1;
end

% --- Insert offsets ---
for idx = 1:length(offset_samples)
    EEG.event(end+1).type     = 'block_offset';
    EEG.event(end).latency    = offset_samples(idx);
    EEG.event(end).urevent    = length(EEG.event) + 1;
end

% --- Sort by latency and ensure consistency ---
[~, sortIdx] = sort([EEG.event.latency]);
EEG.event = EEG.event(sortIdx);
EEG = eeg_checkset(EEG, 'eventconsistency');

end

EEG = insert_block_events(EEG, onset_samples, offset_samples);

%% Actually remove

EEG_clean = EEG;
EEG_clean = pop_select(EEG_clean, 'notime', all_remove);

EEG = eeg_checkset(EEG_clean);

%% Make sure trial events were correctly loaded to EEG structure 

EEG.event.type

for i = 1:min(100, length(EEG.event))
    fprintf('Event %d: type=%s, latency=%.1f\n', ...
        i, EEG.event(i).type, EEG.event(i).latency);
end

% Inspect only block_onset / block_offset events
match_count = 0;

for i = 1:length(EEG.event)
    etype = EEG.event(i).type;

    if ischar(etype) || isstring(etype)
        if strcmpi(etype, 'block_onset') || strcmpi(etype, 'block_offset')
            match_count = match_count + 1;
            fprintf('Block event %d: type=%s, latency=%.1f\n', ...
                match_count, etype, EEG.event(i).latency);
        end
    end
end


%% Remove bad channels 
%specifically looking for channels that do not follow pattern of power
%spectra pl
figure; pop_spectopo(EEG, 1, [], 'EEG' , 'percent',15,'freq',...
        [10 20 30],'freqrange',[1 70],'electrodes','off');

% Signal
%pop_eegplot( EEG_clean, 1, 1, 1);

disp('----------------------')
disp('----------------------')
disp('Manual Input Required')
disp('----------------------')
disp('----------------------')
bad = input('Select bad channels:');

EEG=pop_select(EEG, 'nochannel',bad);

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
EEG = pop_runica(EEG,'icatype','picard','pca',50);

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
% UPDATE file path and naming

EEG = pop_saveset(EEG, 'filename', 'TReC_12454_postEEG_EI.set', 'filepath', ['/Users/elizabethkaplan/Desktop/Greene/subject_folders/TReC_12454']);
 