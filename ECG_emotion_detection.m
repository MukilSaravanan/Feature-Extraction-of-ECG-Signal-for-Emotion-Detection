tStart=tic;

fprintf("Loading the dataset...\n");
if ~exist('DREAMER','var')
    load("C:\Users\user\Documents\MATLAB\DREAMER.mat");
end
fprintf("Dataset got loaded\n");

stimuli_features=[];
baseline_features=[];

Fs=DREAMER.ECG_SamplingRate;
samples_per_min=Fs*60;

total_person=length(DREAMER.Data);
total_stimuli=length(DREAMER.Data{1, 1}.ECG.stimuli);
total_baseline=length(DREAMER.Data{1, 1}.ECG.baseline);
total_channel=size(DREAMER.Data{1, 1}.ECG.stimuli{1,1},2);

for person=1:total_person
    tStartPerson=tic;
    fprintf("Extracting Features from person %d ...\n",person);
    fprintf("Extracted from stimulus ");

for stimulus=1:total_stimuli
    for channel=1:total_channel
        temp_obj=Feature_extracter(DREAMER.Data{1,person}.ECG.stimuli{stimulus,1}(end-samples_per_min+1:end,channel),Fs);
        temp_obj.init();
        stimuli_features(((person-1)*total_stimuli)+stimulus,:,channel)=temp_obj.get_features();
    end
    fprintf("%d \t",stimulus);
end


fprintf("\nExtracted from baseline ");
for baseline=1:total_baseline
    for channel=1:total_channel
        temp_obj=Feature_extracter(DREAMER.Data{1,person}.ECG.baseline{baseline,1}(end-samples_per_min+1:end,channel),Fs);
        temp_obj.init();
        baseline_features(((person-1)*total_baseline)+baseline,:,channel)=temp_obj.get_features();
    end
    fprintf("%d \t",baseline);
end

fprintf("\nExtracted from person %d\n",person);
toc(tStartPerson);
fprintf("===============================================================================================\n");
end
fprintf('Total Features extracted from %d persons\n',person);
toc(tStart);
