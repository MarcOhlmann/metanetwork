
## Resubmission
This is a resubmission. After the package check by the CRAN's team, I:
​
* Added the field \value to the Rd files when it was missing
* Wrote TRUE and FALSE instead of T and F

I did not add references since the paper is still on review. I will update as soon as I get the reference.

Best Regards,

Marc Ohlmann

## Review


Thanks,

If there are references describing the methods in your package, please add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for auto-linking. (If you want to add a title as well please put it in quotes: "Title")

Please add \value to .Rd files regarding exported methods and explain the functions results in the documentation. Please write about the structure of the output (class) and also what the output means. (If a function does not return a value, please document that too, e.g. \value{No return value, called for side effects} or similar)
Missing Rd-tags:
     pipe.Rd:  \value

Please write TRUE and FALSE instead of T and F. Please don't use "T" or "F" as vector names.
'T' and 'F' instead of TRUE and FALSE:

Please fix and resubmit.

Best,
Benjamin Altmann 




## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
