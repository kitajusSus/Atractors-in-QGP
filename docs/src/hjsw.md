# TUTAJ KLIKAĆ DO HJSW



```julia 
dane = load_dataset("datasets/hjsw_lpca/HJSW_DATASET_MAIN.h5")

    tau   = dane[:, 1]
    T_MeV = dane[:, 2]
    A     = dane[:, 3]
    B     = dane[:, 4]

    # Wyliczenie w(t) = tau * T_fm
    w = tau .* (T_MeV ./ MEV_PER_FM)

    # 1. Tworzymy nowy dataset z kolumnami (w, A, B)
    dane_wAB = hcat(w, A, B)

    # 2. Lub dataset z zachowaniem czasu (tau, w, A, B)
    dane_tau_wAB = hcat(tau, w, A, B)


```



i tak o się robi potem robie np. 



```julia 


plot_local_lpca(dane_tau_wAB ).... # i tak dalej 
```
