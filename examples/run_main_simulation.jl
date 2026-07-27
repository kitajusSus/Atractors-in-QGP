using AttractorsQGP

function main()
    model = HJSWModel()

    output_path = "datasets/hjsw_phase_space_1000.hd5"

    println("Uruchamianie symulacji przez run_main...")
    println("Liczba ewolucji (trajektorii): 1000")
    println("Zakres T_0: 100 MeV ... 1000 MeV")
    println("Zakres A_0: -1 ... 6")
    println("Zakres B_0: -1 ... 1")
    println("Jednostki: fermi=1 (T przeliczone z MeV na 1/fm)")

    res = run_main(
        model;
        n_points = 1000,
        tspan = (0.2, 3.0),
        T_range = (100.0, 1000.0),
        A_range = (-1.0, 6.0),
        B_range = (-1.0, 1.0),
        temperature_unit = :MeV,
        saveat = 0.01,
        seed = 42,
        output_file = output_path
    )

    println("Symulacja zakończona!")
    println("Kształt wygenerowanej macierzy przestrzeni fazowej [tau, T, A, B]: ", size(res.dataset))
    println("Plik zapisano w: ", output_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
