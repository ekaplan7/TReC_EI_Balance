% code adapted from existing preproc pipeline from Dr. Danielle Sliva
clear all;
close all;
clc;

% eeglab
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

%% SETUP

subjectID = '12332';     % change per subject
session   = 'post';      % 'pre' or 'post'

%rawdatadir   = '/Users/elizabethkaplan/Desktop/Greene/source_data';
subject_dir = fullfile('/Users/elizabethkaplan/Desktop/Greene/subject_folders',['TReC_' subjectID]);
analysis_dir = '/Users/elizabethkaplan/Desktop/Greene/TReC_EI_Balance';

% Main file to preprocess
setfile = ['TReC_' subjectID '_' upper(session(1)) session(2:end) '.set'];

% Output filename
output_setfile = ['TReC_' subjectID '_' upper(session(1)) session(2:end) 'EEG_EI_NEW.set'];

%% LOAD DATA

EEG = pop_loadset('filename', setfile, 'filepath', subject_dir);

%seelct only eeg channels 
EEG = pop_select(EEG, 'channel', find(strcmpi({EEG.chanlocs.type}, 'eeg')));

% drop non-EEG channels
EEG = pop_select(EEG, 'nochannel', find(ismember(lower({EEG.chanlocs.labels}), {'eog','ekg'})));

% drop any channel whose label starts with BIP
EEG = pop_select(EEG, 'nochannel', find(startsWith({EEG.chanlocs.labels}, 'BIP', 'IgnoreCase', true)));

% Save in ALLEEG for later visual inspection
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 1);

chanloc_file = '/Users/elizabethkaplan/Desktop/eeglab2026.0.0/sample_locs/standard_waveguard64.elc';
EEG = pop_chanedit(EEG, 'load', {chanloc_file, 'filetype', 'autodetect'});
EEG_mask = EEG; 

elc_locs = readlocs(chanloc_file);

n_matched = 0;
for iChan = 1:length(EEG.chanlocs)
    idx = find(strcmpi(EEG.chanlocs(iChan).labels, {elc_locs.labels}));
    if ~isempty(idx) && ~isempty(elc_locs(idx(1)).X) && ~isnan(elc_locs(idx(1)).X)
        EEG.chanlocs(iChan).X = elc_locs(idx(1)).X;
        EEG.chanlocs(iChan).Y = elc_locs(idx(1)).Y;
        EEG.chanlocs(iChan).Z = elc_locs(idx(1)).Z;
        n_matched = n_matched + 1;
    else
        fprintf('No coordinate match for: %s\n', EEG.chanlocs(iChan).labels);
    end
end
fprintf('Matched %d of %d channels\n', n_matched, length(EEG.chanlocs));

EEG.chanlocs = convertlocs(EEG.chanlocs, 'cart2all');  % fills in theta/radius/sph from X/Y/Z
EEG = eeg_checkset(EEG);

%% ADD IN BLOCKS (REST VS SUPPRESS)

block_file = fullfile(subject_dir, ...
    ['study-353_sub-' subjectID '_task-eeg_ses-' session '.csv']);

block_table = readtable(block_file);

% Use columns 2 and 3 (safer than relying on names)
Blocks_onset_ms  = block_table{:,2}';
Blocks_offset_ms = block_table{:,3}';

% sanity check
if length(Blocks_onset_ms) ~= length(Blocks_offset_ms)
    error('Block onset and offset columns mismatch');
end

fs = EEG.srate;
video_delay_sec = 1.26690;

% Convert ms -> sec
block_onsets_sec  = Blocks_onset_ms  / 1000;
block_offsets_sec = Blocks_offset_ms / 1000;

% Add video delay
block_onsets_sec  = block_onsets_sec  + video_delay_sec;
block_offsets_sec = block_offsets_sec + video_delay_sec;

% Convert sec -> samples
onset_samples  = block_onsets_sec  * fs;
offset_samples = block_offsets_sec * fs;

% Force row vectors
onset_samples  = onset_samples(:)';
offset_samples = offset_samples(:)';

%% PREPROCESS DATA 

%remove DC offset 
for iChan = 1:size(EEG.data, 1)
    EEG.data(iChan, :) = single(double(EEG.data(iChan, :)) - ...
                                mean(double(EEG.data(iChan, :))));
end

%Notch filter 
EEG_params.notch.low_cutoff  = 58;
EEG_params.notch.high_cutoff = 62;
notch_filter = 1;

EEG = pop_eegfiltnew(EEG, EEG_params.notch.low_cutoff, ...
    EEG_params.notch.high_cutoff, [], notch_filter);
EEG = eeg_checkset(EEG);

% high & low pass filter
EEG = pop_basicfilter(EEG, 1:EEG.nbchan, ...
    'Cutoff', 0.1, ...
    'Design', 'butter', ...
    'Filter', 'highpass', ...
    'Order', 2, ...
    'RemoveDC', 'on', ...
    'Boundary', 'none');

EEG = pop_basicfilter(EEG, 1:EEG.nbchan, ...
    'Cutoff', 70, ...
    'Design', 'butter', ...
    'Filter', 'lowpass', ...
    'Order', 4, ...
    'Boundary', 'none');

%Run amplitude sliding window 
fs = EEG.srate;
winSize_sec = 0.05;
step_sec    = 0.025;
threshold   = 200;
channels    = 1:EEG.nbchan;

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

% padding around detected artifact windows
pad_sec = 0.05;
padSamp = round(pad_sec * fs);

if padSamp > 0
    rejMask = conv(double(rejMask), ones(1, 2*padSamp + 1), 'same') > 0;
end

% Convert mask to start/end sample ranges
diffMask = diff([0, rejMask, 0]);
starts = find(diffMask == 1);
ends   = find(diffMask == -1) - 1;

% Convert to seconds
artifactRanges_sec = [starts(:) ends(:)] / fs;

fprintf('Built %d sliding-window artifact periods\n', size(artifactRanges_sec,1));

%% MERGE ARTIFACT RANGES (tics and artifact) 

all_remove = artifactRanges_sec;

if isempty(all_remove)
    warning('No artifact ranges to remove.');
    all_remove = zeros(0,2);
else
    all_remove = sortrows(all_remove, 1);
    merged = all_remove(1,:);
    for i = 2:size(all_remove,1)
        if all_remove(i,1) <= merged(end,2)
            merged(end,2) = max(merged(end,2), all_remove(i,2));
        else
            merged = [merged; all_remove(i,:)]; 
        end
    end
    all_remove = merged;
end

%% PROTECT BLOCK ONSET AND OFFSET MARKERS 

protect_margin = 0.004; % 4 ms margin
original_remove = all_remove;

for i = 1:length(block_onsets_sec)
    for j = 1:size(original_remove, 1)

        % Protect onset
        if block_onsets_sec(i) >= original_remove(j,1) && block_onsets_sec(i) <= original_remove(j,2)
            fprintf('Onset %.3f s inside removed segment (%.3f-%.3f s). Trimming.\n', ...
                block_onsets_sec(i), original_remove(j,1), original_remove(j,2));
            all_remove(j,2) = block_onsets_sec(i) - protect_margin;
        end

        % Protect offset
        if block_offsets_sec(i) >= original_remove(j,1) && block_offsets_sec(i) <= original_remove(j,2)
            fprintf('Offset %.3f s inside removed segment (%.3f-%.3f s). Trimming.\n', ...
                block_offsets_sec(i), original_remove(j,1), original_remove(j,2));
            all_remove(j,1) = block_offsets_sec(i) + protect_margin;
        end
    end
end

% Remove invalid windows created by trimming
if ~isempty(all_remove)
    all_remove = all_remove(all_remove(:,2) > all_remove(:,1), :);
end

% Re-merge after trimming
if ~isempty(all_remove)
    all_remove = sortrows(all_remove, 1);
    merged = all_remove(1,:);
    for i = 2:size(all_remove,1)
        if all_remove(i,1) <= merged(end,2)
            merged(end,2) = max(merged(end,2), all_remove(i,2));
        else
            merged = [merged; all_remove(i,:)]; 
        end
    end
    all_remove = merged;
end

fprintf('Final number of protected artifact periods to remove: %d\n', size(all_remove,1));

%% INSERT BLOCK EVENTS

EEG = insert_block_events(EEG, onset_samples, offset_samples);
EEGfull = EEG;

%% REMOVE ARTIFACTS

EEG_clean = EEG;

if ~isempty(all_remove)
    EEG_clean = pop_select(EEG_clean, 'notime', all_remove);
else
    warning('No segments removed.');
end

EEG = eeg_checkset(EEG_clean);

%% CHECK BLOCK EVENT MARKERS - make sure all expected markers are present in data 

disp('EEG.event.type values:');
for i = 1:length(EEG.event)
    disp(EEG.event(i).type);
end

for i = 1:min(100, length(EEG.event))
    fprintf('Event %d: type=%s, latency=%.1f\n', ...
        i, string(EEG.event(i).type), EEG.event(i).latency);
end

match_count = 0;
for i = 1:length(EEG.event)
    etype = EEG.event(i).type;

    if ischar(etype) || isstring(etype)
        if strcmpi(etype, 'block_onset') || strcmpi(etype, 'block_offset')
            match_count = match_count + 1;
            fprintf('Block event %d: type=%s, latency=%.1f\n', ...
                match_count, string(etype), EEG.event(i).latency);
        end
    end
end

%% REMOVE BAD CHANNELS

%remind EEGLAB of chanlocs (idk why but will through a topo error if not
%ran again)
chanloc_file = '/Users/elizabethkaplan/Desktop/eeglab2026.0.0/sample_locs/standard_waveguard64.elc';

EEG = pop_chanedit(EEG, 'lookup', chanloc_file);
EEG = eeg_checkset(EEG);

% Also try converting location formats
EEG.chanlocs = convertlocs(EEG.chanlocs, 'cart2all');

figure;
pop_spectopo(EEG, 1, [], 'EEG', 'percent', 15, 'freq', [10 20 30], ...
    'freqrange', [1 70], 'electrodes', 'off');

disp('----------------------')
disp('Manual Input Required')
disp('----------------------')
bad = input('Select bad channels (indices in CURRENT EEG): ');

bad = bad(:)'; % force row vector
bad = bad(bad >= 1 & bad <= EEG.nbchan); % keep only valid

EEG = pop_select(EEG, 'nochannel', bad);

%% RE-REFERENCE
EEG = pop_reref(EEG, []);  % average reference (excludes already-removed bad channels)
EEG = eeg_checkset(EEG);

%% ICA

% High-pass at 1 Hz before ICA
EEG = pop_basicfilter(EEG, 1:EEG.nbchan, ...
    'Cutoff', 1, ...
    'Filter', 'highpass', ...
    'Design', 'butter', ...
    'Order', 2);

% Downsample first for faster ICA
EEG = pop_resample(EEG, 256);

% Run ICA
% Make sure PCA number is < number of channels after bad-channel removal
ica_pca_n = min(50, EEG.nbchan - 1);
EEG = pop_runica(EEG, 'icatype', 'picard', 'pca', ica_pca_n);

%% VISUALIZE AND REJECT IC COMPONENTS

figure;
pop_selectcomps(EEG, 1:min(35, size(EEG.icaweights,1)));
pop_eegplot(EEG, 0, 1, 1);

disp('----------------------')
disp('Manual Input Required: Select Components to Reject')
disp('----------------------')

reject_components = input('Enter the list of components to reject (e.g., [1 2 5 8]): ');

EEG = pop_subcomp(EEG, reject_components, 0);

%% INTERPOLATION

EEG = pop_interp(EEG, EEGfull.chanlocs, 'spherical');

%% FINAL VISUAL INSPECTION

EEG.icaweights = [];
EEG.icasphere = [];
EEG.icawinv = [];
EEG.icaact = [];
EEG.icachansind = [];

EEG = eeg_checkset(EEG);
pop_eegplot(EEG, 1, 1, 1);
eegh;

%% SAVE DATA 

output_setfile = ['TReC_' subjectID '_' upper(session(1)) session(2:end) 'EEG_EI.set'];

EEG = pop_saveset(EEG, 'filename', output_setfile, 'filepath', subject_dir);

%% 
%% ------------------------------------------------------------------------
% LOCAL FUNCTION
% -------------------------------------------------------------------------

function EEG = insert_block_events(EEG, onset_samples, offset_samples)
% insert_block_events Add trial onset/offset events to EEG struct
%
% Inputs:
%   EEG             - EEGLAB dataset struct
%   onset_samples   - vector of block onset latencies (in samples)
%   offset_samples  - vector of block offset latencies (in samples)
%
% Output:
%   EEG             - EEGLAB dataset struct with added events

if length(onset_samples) ~= length(offset_samples)
    error('Onset and offset vectors must be the same length.');
end

% Initialize EEG.event if empty
if ~isfield(EEG, 'event') || isempty(EEG.event)
    EEG.event = struct('type', {}, 'latency', {}, 'urevent', {});
end

% Insert onsets
for idx = 1:length(onset_samples)
    EEG.event(end+1).type    = 'block_onset';
    EEG.event(end).latency   = onset_samples(idx);
    EEG.event(end).urevent   = length(EEG.event);
end

% Insert offsets
for idx = 1:length(offset_samples)
    EEG.event(end+1).type    = 'block_offset';
    EEG.event(end).latency   = offset_samples(idx);
    EEG.event(end).urevent   = length(EEG.event);
end

% Sort by latency
[~, sortIdx] = sort([EEG.event.latency]);
EEG.event = EEG.event(sortIdx);

EEG = eeg_checkset(EEG, 'eventconsistency');
end