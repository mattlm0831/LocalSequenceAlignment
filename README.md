# LocalSequenceAlignment
Github Repo with code from a research project I completed that explores adapting the Smith-Waterman DP algorithm to GPU architecture.


This program occurred in several sequences, namely: sequential C, Multi-threaded C, and the CUDA C with hyper-parallelism. It was 
done in this way to figure out the original concurrency road blocks associated with parallelism before venturing into the uncharted
GPU territory. As you will find in this repo, many interesting concepts were employed to overcome these road blocks.

We were able to attain considerable speed-ups from the multi-threaded version, especially when considering data transfer rates.

Samples were actual sequences, smaller was about 100 proteins in length, and the larger was an entire chromosome sample.

Please let me know what you think.
