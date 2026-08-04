# NCBJ docs

Bardzo przepraszam ale w ramach zwiększenia zrozumienia kodu i ideii nazwy są pisane
dwunastozgłoskowcem

## Jak używać `example/ncbj_lle.jl` ?
> wszystko przygotowane pod prace z repl

- krok 1
uruchamianie środowiska
```julia
using AttractorsQGP
includet("examples/ncbj_lle.jl")
```

- krok 2 Tworzenie danych do  testowania

```julia
wektor_x = [0.2:5...]
 macierz = ncbj1_macierz_wszyskich_punktów!(wektor_x)
```


- krok 3 znajdowanie sąsiadów dla wybranego punktu
```julia
ile_sąsiadów = 3
indeks_wybranego_punktu = 3
sasiedzi, indeksy_sasiadow, wybrany_punkt = ncbj2_sąsiedzi!(macierz, indeks_wybranego_punktu, ile_sąsiadów )
```

> należy pamietać że nazwy służą  temu by każdy wiedział co jest czym bez szczególnego zastanawiania


- krok 4  obliczanie macierzy wag


```julia
ncbj3_calculate_wagi_dla_x_i!(sasiedzi, wybrany_punkt, ile_sąsiadów)




```

## Dodane nowe funkcje i porządna "refactoryzacja"

jak uruchamiać przykłady? 

```julia
julia src/examples/lle_examples/swissroll.jl
julia src/examples/lle_examples/helix.jl
```


# `examples/local_pca/lpca.jl`




```julia
includet("examples/local_pca/ex_lpca.jl")
dane = load_dataset("datasets/....")


ex_lpca(dane2, 10,10)
Tau | Śr. wymiar dla K=10 | Śr. wymiar dla K=20 | Zmiana (%)
---------------------------------------------------------------------------
0.2 | 2.0 | 2.0 | 0.0%
1.18 | 1.022 | 1.0 | -2.152641878669278%
2.16 | 1.0 | 1.0 | 0.0%
3.14 | 1.0 | 1.0 | 0.0%
4.12 | 1.0 | 1.0 | 0.0%
5.1 | 1.0 | 1.0 | 0.0%
6.08 | 1.0 | 1.0 | 0.0%
7.0600000000000005 | 1.0 | 1.0 | 0.0%
8.040000000000001 | 1.0 | 1.0 | 0.0%
9.020000000000001 | 1.0 | 1.0 | 0.0%

 ex_lpca(dane, 10,10)

Tau | Śr. wymiar dla K=10 | Śr. wymiar dla K=20 | Zmiana (%)
---------------------------------------------------------------------------
0.22 | 2.0 | 2.0 | 0.0%
0.31 | 2.0 | 2.0 | 0.0%
0.4 | 1.902 | 1.821 | -4.258675078864351%
0.49 | 1.572 | 1.511 | -3.8804071246819443%
0.58 | 1.362 | 1.319 | -3.1571218795888507%
0.67 | 1.216 | 1.16 | -4.605263157894741%
0.76 | 1.113 | 1.061 | -4.672057502246186%
0.85 | 1.039 | 1.0 | -3.7536092396535063%
0.94 | 1.004 | 1.0 | -0.39840637450199234%
1.03 | 1.0 | 1.0 | 0.0%

```




