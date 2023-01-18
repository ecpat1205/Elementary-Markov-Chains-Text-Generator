# Elementary-Markov-Chains-Text-Generator

## Markov Chain Text Generator Design
### Transition Matrix Design
In the generator, the characters of an original text are simplified, followed by the generated one. In the original text, letters are turned to lower case, and punctuation marks are replaced with the symbol "\_", so for a first-order transition matrix, the states of a Markov Chain are the 26 English lower-case letters and the symbol "\_". Clearly, the first-order transition matrix is a square matrix, but note that here the matrix with transition order greater than 1 is non-standard since it is not square, while every entry represents the probability with every row add up to one. This situation is based on the design of the transition probability matrix. In the matrix, each row index is the existing k-gram state in the original text, and each column index is one of the 26 lower-case letters and "\_". For example, in a tiny text "Hello world", the 3-gram states are "hel", "ell", "llo", "lo ", "o w", " wo", "wor", "orl" and "rld". In the original text from which we want to simulate, existing states are extracted instead of doing the permutation for all letters and symbols to even generate states never appearing, which tries to save the memory when running the generator program. The following matrix is an example of one third-order transition matrix in this generator constructed.

Two problems are met when generating characters from the transition probability matrix designed here. The first problem is that the probability of a existing state moving to any other character may be non-existent in the original text. If the state only exists at the end of the text, then the relative frequency would be NaN. Generally it always happens when the transition order is greater than 1 if it happens when order is 2. It could happen as well when order is 1, e.g. the only "." is at the end of the text when it is not replaced with "_". This row of the matrix is detected and removed, because it is equivalent to the non-existent state
with respect to generating a character, but it is considered in the initial state probability distribution still in case needed for initial states generation. The second problem is dealing with the character generation from never-occurring states in the original text that are generated by the last state, to which the following section gives the solution.

### Text Generator Design
For a simple text generator, when the current k-gram state does not exist in the row index of the transition matrix, the next character is generated simply based on the single-character probability distribution. For a more complete text generator, the next character is expected to be generated by its previous $(k-1)$ characters. If the new state does not exist in the $(k-1)$-order transition matrix as well, then the previous $(k-2)$ characters are used. Recursively, the next character is generated until the previous $(k-m)$ characters state exist in the $(k-m)$-order transition matrix where $m$ refers to the order reduction times. If unfortunately when $k=1$ the state does not appear in the original text, then the next character would be generated by the single-character probability distribution.

## Experiment
### Experiment Design
The original text used for simulation is the Inaugural Address of Barack Obama as the president in the US in 2009, which can be found via http://obamaspeeches.com/P-Obama-Inaugural-Speech-Inauguration.htm. One-tenth, half, and full of the characters in the original text, which are 1337, 6687, 13374 characters respectively,  are put into the simple and less simple generators with three transition orders $1$, $4$, $8$ to generate a 1200-grams long text respectively, with 26 lower-case letters and the symbol "\_". Under each condition, the simulation runs 30 times and collects the word generated accuracy, the average time per character generation, the initial transition matrix making time, the average total time of generating a text. The word accuracy is calculated by total counts of real words divided by total counts of letters combination generated. The symbol "\_" is removed when summarizing the counts. For the complete generator, the average of loops caused by the transition order reduction when generating a character is also included.
### Demonstration of Some Results
The outputs of simulations of two generators by setting $4$ as the transition order with one-tenth and full of the original text characters are displayed here.
#### Simple Generator with 10\% Characters
The text generated is displayed in Figure 1. The plot of generation time per characters is shown in Figure 2. The total processing and generating time are 5.683 seconds. The failure to generate real words at the end of the text coincides with the jump of generation time when generating a character later than 800th, which means the generator meets with the situation when some states are not in the transition matrix, and it has to randomly produce some characters based on the single-character distribution.

![Alt text](/Figures/short_simple_four_text.png "Figure 1: Text (Simple generator with 10\% characters)")