# Elementary-Markov-Chains-Text-Generator

## Markov Chain Text Generator Design
### Transition Matrix Design
In the generator, the characters of an original text are simplified, followed by the generated one. In the original text, letters are turned to lower case, and punctuation marks are replaced with the symbol "\_", so for a first-order transition matrix, the states of a Markov Chain are the 26 English lower-case letters and the symbol "\_". Clearly, the first-order transition matrix is a square matrix, but note that here the matrix with transition order greater than 1 is non-standard since it is not square, while every entry represents the probability with every row add up to one. This situation is based on the design of the transition probability matrix. In the matrix, each row index is the existing k-gram state in the original text, and each column index is one of the 26 lower-case letters and "\_". For example, in a tiny text "Hello world", the 3-gram states are "hel", "ell", "llo", "lo ", "o w", " wo", "wor", "orl" and "rld". In the original text from which we want to simulate, existing states are extracted instead of doing the permutation for all letters and symbols to even generate states never appearing, which tries to save the memory when running the generator program. The following matrix is an example of one third-order transition matrix in this generator constructed.