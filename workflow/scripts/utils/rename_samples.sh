

rename 's/NL4/R/; s/NL1/F/' *.ab1

rename 's/ (F|R)_[A-Z0-9]+/_$1/' *.ab1

# if necessary:
rename 's/^M/new_/' *.ab1