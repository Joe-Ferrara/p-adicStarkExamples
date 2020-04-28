# p-adic Stark Examples

This is the code used to compute the examples in the paper, [<i>A p-adic Stark conjecture in the rank one setting</i>](https://arxiv.org/abs/1904.10561), and my thesis, [<i>Stark's Conjectures for p-adic L-functions</i>](https://escholarship.org/uc/item/4qv5b8tz#main).

The folders ExQ17p5Prec60, ExQneg23p5Prec60, and ExQneg31p3Prec60, contain code written to compute the examples for the conjecture at s = 0 and s = 1. The folders ExQ17p5Prec63, ExQneg23p5Prec72, ExQneg31p3Prec77 contain code that reproduces and reaffirms the calculations in ExQ17p5Prec60, ExQneg23p5Prec60, and ExQneg31p3Prec60 for the conjecture at s = 0.

The folder OverconvergentModularSymbols contains code written by Rob Pollack and Rob Harron to compute overconvergent modular symbols (OMS). The code in the examples folders uses the code in the overconvergent modular symbols folders to compute values of a p-adic L-function. More extensive OMS code from which the OMS code here was taken from can be found at Rob Pollack's github page: https://github.com/rpollack9974/OMS. The only significant new code in the OMS folder here is in the pLfunctionNew.sage file and the Functions.sage file.
