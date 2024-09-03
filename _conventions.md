# Math notation

- Decimals: write $0.2$, not $.2$ (non-standard AFAIK--SV comment).
- All distributions are written \mathit{Distribution} in sampling statements
- All distributions are written as distribution in the text except for distributions that are proper nouns such as Poisson and Bernoulli.
- log-normal:
    - log-normal in the text
    - \mathit{LogNormal} in sampling statements 
- \Theta for arbitrary parameters
- \theta for bernoulli/binomial
- \boldsymbol{\theta} for categorical, multinomial
- y for arbitrary data
- all vectors and parameters are in bold
- elements of a vector should be inside [...,...]
- iterate over observations with n = {1,2,..,N}
- iterate over subjects with i = {1,2,..,N_subj}
- iterate over items with j = {1,2,..,N_item}
- we use the Gelman/stan notation: u[subj[nrow]] to indicate the adjustment that corresponds to the subject subj[nrow]
- we start from 1 rather than 0.
- always write i=1,...,N not i=1...N
- always write \exp and \log, never exp and log, in math mode
- we use $4000$ rather than $4,000$

# Code:

Follows this https://style.tidyverse.org/

* avoid using `round` in the output, Rmarkdown already does this for us

# Variable conventions for the code

We need to work on something like this (This can change a lot.):

* dataframe with data of Grodner and Gibson 2008

`df_gg08`

* dataframe with fake data

`df_fake`


* list with data of Grodner and Gibson 2008

`lst_gg08`

* Stan fit with a lognormal model of data of Grodner and Gibson 2008

`fit_lognormal_gg08`

* dataframe version of the Stan fit with a lognormal model of data of Grodner and Gibson 2008

`df_fit_lognormal_gg08`

* list version of the Stan fit with a lognormal model of data of Grodner and Gibson 2008

`lst_fit_lognormal_gg08`

# Linguistic examples

Need to decide and be consistent on the following:
* italics?
* use of "
* both?
* what about a linguistic example in, say, German, and its translation?

# Conventions for text

- We refer to functions as ``` `functionname()` ```

# English

* Do not use "we" so much, especially not in code comments, but also not in text.
* Always say "replace X with Y", not "replace X for Y".
* Space after ms: not 500ms, but 500 ms.
* Pay attention to the difference between its and it's.
* E.g., and i.e., etc., have to be followed by a comma.
* Figure is always Figure, not Fig, Fig., etc. (I have also made this mistake, am fixing it now).
* Never use "We here". This is not idiomatic English.
* Never use "allows to". This is not English. It has to be "allows us to" or "allows the researcher to", etc.
* Never use "Note that" unless there is really no other way to say this.
* "to the contrary" is not really used correctly, and it is also not idiomatic English. Just by "By contrast", and always begin a sentence with "By contrast", don't put it in the middle of the sentence; placing it initially prepares the reader for a contrast and improves readability.
* data set rather than dataset or data-set
* data frame rather than data-frame, and dataframe
* cross-validation rather than cross validation.
* response times rather than reaction times
* can not does not mean the same thing as cannot; we use "can not" incorrectly in the book. (I fixed these cases).
* Use learned and not learnt.
* varying intercepts (plural) and varying slopes (plural) model

* Punctuation should always be contained within single or double quotes, e.g., when ending a sentence. So, write: "this is the end of a sentence." and not "this is the end of a sentence".
* Do not italicize whole phrases or sentences for emphasis; use syntactic reformulation to emphasize. Example:

Bad: *It is still not the appropriate model for the generative process that we are assuming in this section*, but it is a good start.

Good: It is still not the appropriate model for the generative process that we are assuming in this section, but it is a good start.

# Parts of the book
We use section, chapter, Figure, Table, Box, and Equation before a reference `\@ref(something)`.


