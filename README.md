# Feature Extraction of ECG Signal for Emotion Detection
## Instructions 
- To extract features from the complete ECG dataset of DREAMER dataset, run 'ECG_emotion_detection.m' in matlab. The extracted stimuli and baseline features will be stored in the respective variables 'stimuli_features' and 'baseline_features' which are saved as 'Extracted_features414x17x2.mat' for training the neural network in the next stage.
- For testing, use 'test_code.m', it works the same but with a smaller DREAMER dataset 'ECG_sample_dataset.mat'
- The live script 'affective_dimensional_model.mlx' gives a walkthrough of the operations

### Note
1) The class definition file 'Feature_extracter.m' should be in the same working directory.
2) The complete dataset can be downloaded from the following link https://zenodo.org/record/546113#.YS0TQFvhWV4
 
## Best ECG Datasets for emotion detection

- AMIGOS dataset: 
www.eecs.qmul.ac.uk/mmv/datasets/amigos/index.html, 
https://github.com/pokang-liu/AMIGOS

- DREAMER dataset:
https://zenodo.org/record/546113#.YS0TQFvhWV4

- AuBT:
https://www.informatik.uni-augsburg.de/de/aktuell/kolloquium/previous/ws07-08/2008-02-07kim/index.html

- ASCERTAIN:
https://ascertain-dataset.github.io/

## References 

- The presentation slide 'AI based Emotion Detection with ECG signal.pdf' gives more information on the flow
- Some of the reference articles can be found in the child directory