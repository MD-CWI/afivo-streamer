# Reaction set

In this file, two reaction sets are stored.

# Komuro's reaction set

A reaction set for N2-O2 mixtures from Ono and Komuro 2020 (https://doi.org/10.1088/1361-6463/ab4e65) is presented in `Komuro_reactions.csv`.

# Baohong's reaction set

We constructed a new reaction set for N2-O2 mixtures in `Baohong_reactions.csv`. This is primarily from:

Kossyi 1992 (https://doi.org/10.1088/0963-0252/1/3/011)

Gordiets 1995 (https://doi.org/10.1109/27.467998) 

Aleksandrov 1999 (https://doi.org/10.1088/0963-0252/8/2/309)

Fresnet 2002 (https://doi.org/10.1088/0963-0252/11/2/305)

Tochikubo 2002 (https://doi.org/10/d9tq56)

Atkinson 2004 (https://doi.org/10.5194/acp-4-1461-2004)

Florescu 2006 (https://doi.org/10.1016/j.physrep.2006.04.002)

Gordillo 2008 (https://doi.org/10.1088/0022-3727/41/23/234016)

Popov 2011 (https://doi.org/10.1088/0022-3727/44/28/285201)

Panchesnyi 2013 (https://doi.org/10.1088/0022-3727/46/15/155201)

# Convert reaction format

The csv format of a reaction set can be converted to the format used in `Afivo-streamer` via `chemistry_reaction_parser.py`:

    python chemistry_reaction_parser.py Baohong_reactions.csv -convert_tex Baohong_converted.csv 
    python chemistry_reaction_parser.py Baohong_converted.csv
