## CRIMAP

Code taken from <http://compgen.rutgers.edu/crimap>

In particular <http://compgen.rutgers.edu/downloads/crimap.source.zip>

Documentation at <https://saf.bio.caltech.edu/saf_manuals/crimap-doc.html>

Compile the software with:

```
gcc -O -o crimap *.c -lm
```

### Changes to code

Made the following changes to get it to compile:

- In `e_ped.c`, added argument to calls to `exit()`

- In `our_allo.c` and `our_orde.c`, moved definition of `morecore` to
  earlier in the file

- In `defs.h`, changed `typedef INT` from `int` to `long`

- Adjusted a bunch of parameters in `defs.h`, `e_conv.c`, and
  `e_crimap.c`

- In `e_crimap.c`, changed print statements to print estimated map
  with greater precision
