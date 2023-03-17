## Resubmission
This is a resubmission. The previous submission was rejected because of a function that
writes in the user's home filespace.
However, the offending function lives under the "tools/" folder and is only called during
the package configuration to write a "Makevars" file in the "src/" folder that selects the appropriate 
C++ version. My understanding is that this should be in line with CRAN policies. I apologize if I 
misunderstood and this is not the case.

## R CMD check results

0 errors \| 0 warnings \| 0 note

This is a new release.
