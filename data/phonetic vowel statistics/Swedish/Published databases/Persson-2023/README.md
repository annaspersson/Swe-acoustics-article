This folder contains txt and csv files with vowel data of four types. 1) Vowel data generated by the Praat vowel analyzer script (/raw-formant-data). 2) Vowel data that has been outlier corrected by the Persson-2021-outlierCorrection.rmd script (/outlier-corrected-vowel-data). 3) Vowel statistics generated by the get_vowel_statistics_Swe.rmd script (/vowel_statistics), with the main purpose of using it for resynthesis. 4) Vowel statistics transformed into different cue spaces, generated by the Persson-2021-repair.rmd script (/vowel_statistics).
The scripts relate to each other in the following way: the outlierCorrection script reads in the raw formant data generated by the Praat script and outputs csv-files that are read in by the get_vowel_statistics script. The get_vowel_statistics script outputs vowel statistics that are read in by the repair script, that adds and repair vowel formant information and outputs csv-files that are read in by thesis chapter rmds and by experiment report rmds.

The subfolder *vowel_statistics* includes the following types of statistics, all generated on the outlier corrected data: general vowel statistics; vowel statistics transformed into different cue spaces by normalization functions; vowel statistics that are used for generating resynthesized stimuli. The latter stats type either allow for random sampling within the praat script, or hand individual tokens to be resynthesized by the praat script. The tokens are generated either by random sampling in R or by generating a test grid through manual inspection.