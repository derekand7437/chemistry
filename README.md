# The Chem Bench

A study site for first-semester (10th grade) chemistry: a 21-day guided path that starts from
"everything is made of atoms" and ends at gas laws and periodic trends, plus a reference library
with practice problems that never run out.

The whole site is one file — `index.html` — with its CSS and JavaScript inside it. No build step,
no dependencies, no server code. Progress saves in each visitor's own browser, so the site
remembers which day they are on.

## The 21 days

Atoms and symbols · metals, nonmetals and ions · ionic vs covalent · naming ionic compounds ·
writing formulas · polyatomic ions · Roman numerals · covalent prefixes · mixed naming ·
counting atoms · balancing equations · combustion · molar mass · grams and moles ·
particles and STP · mole ratios · grams to grams · Boyle and Charles · the ideal gas law ·
periodic trends.

## Publishing it

Push this folder to a GitHub repo, then **Settings → Pages** → *Deploy from a branch*,
branch `main`, folder `/ (root)`. The site appears at
`https://YOUR-USERNAME.github.io/REPO-NAME/`.

## Editing it

Open `index.html`. The daily lessons are the `PLAN` array near the bottom of the `<script>`;
the practice problems come from the generator functions just above it.
