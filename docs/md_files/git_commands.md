# Main git commands

Here you can find essential git commands which you need to work with a github repository.

`git clone "URL"` - make a copy of the repository on your local machine <br>

## Branches

`git branch` - shows the current branch<br>
`git checkout -b BranchName` - creates a new branch<br>
`git branch -d BranchName` - removes branch locally<br>
`git push origin --delete BranchName` - removes branch online<br>

## Changes in your branch
`git status` - shows which files have been changed<br>
`git add` - addes file into repository<br>
`git rm -r Filename` - removes file from repo<br>
`git commit -m "message"` - commits changes<br>
`git push` - pushes your local commits to the server<br>
`git pull` - pulls all updates from server<br>

`git stash save --keep-index && git stash drop` - discards changes which you didn't commit
`git checkout "commit_tag" src/Lagrangian/Lagrangian.c` - replaces file (src/Lagrangian/Lagrangian.c) with the same (file) from previous commit

## Rename branch
Change into branch which should be renamed:<br>
`git checkout old_name`
Rename local branch:<br>
`git branch -m new_name`
Push the new_name local branch and reset the upstream branch:<br>
`git push origin -u new_branch`
Delete the old_name remote branch:<br>
`git push origin --delete old_name` 

