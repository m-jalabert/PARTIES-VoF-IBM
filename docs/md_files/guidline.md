# Workflow and Branching
The workflow is orginized in the following way:
There is a  `master` branch - this is the main branch with features which have been developed and tested. This branch is considered as the most stable
 branch and recomended to be used for simulations.
 
## Code development
Each contributor of the repository is able to develop a new functionality by creating a new branch. It is preferabe that this branch originates 
from the latest version of the master branch.

The new branch also should be named with the following way:<br> 
Last name of a creator + underscore + short flag
which indicates the purpose of the branch, for example `metelkin_develop` or `metelkin_postproc`. 

## Merging your branch with the master branch
In order to minimize the amount of bugs in the master branch the following steps have to be made to merge a new branch with the master branch:

 1. **Pass testcases**. <br> When the implemented code from a new branch is finished, it has to be tested using the collection of test cases to be sure that
the code functionality is not affected. The guideline of how to perform test cases is located here(link).
 2. **Create new testcase** <br> If your code includes an additional simulation tools (for example, new solver) it is strongly recomended to create your own test case
and add it to the list of test cases. The guideline of how to add new test cases is located here (link and the detailed guidleine of how to create a test case).
 3. **Document your changes**. <br> The last step before merging your branch with the master branch is to document the changes which you have implemented in your new branch. The description should be added
in the corresponding chapter of documentation.

If all the steps are completed the new branch should be merged with the master branch via **Pull Request**.

