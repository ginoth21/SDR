# Work Breakdown Disclosure Document

## Group 09

| Member Name | MacID | Student Number |
| --- | --- | --- |
| Ginoth Kumarakulasingam | kumarakg | 400258323 |
| Ty Greenwood | greent14 | 400311174 |
| Oscar | lus70 | 400252202 |
| Cass Smith | smithc85 | 400255184 |

## Contributions

| Member | Mono | Stereo | RDS | Threading |
| --- | --- | --- | --- | --- |
| Ginoth | Block convolution, low-pass filter generator | Worked with Cass on stereo path implementation| Root-raised cosine filter, debugging | Worked on threading for mono/stereo and debugged RDS threading |
| Ty | Python modeling of convolve/downsample and resample, C++ implementation of convolve/downsample and resample (minus state saving), testing | Python modeling of all pass filter, testing | zero/peak detector, Manchester decoder, differential decoder, frame synchronization/traversal and syndrome calculator | Second vector queue for RDS thread |
| Oscar | Convolve/downsample function, fast resampler. same as Ty |same as Ty | Same as Ty | Same as Ty |
| Cass | Stdin/stdout functions, interleaving functions | PLL port, bandpass filter generator, mixer, fast all-pass filter, combiner | IQ PLL, initial signal processing path, PI and PT extraction functions | Vector queue with mutex locks |
