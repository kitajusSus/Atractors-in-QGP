### Źródło Różnic: Różne Teorie Hydrodynamiczne

Oba fragmenty pokazują rozwinięcie anizotropii ciśnień $A(w)$ dla dużych $w$ (rozwinięcie gradientowe). Wyglądają podobnie, ale różnią się współczynnikiem przy członie drugiego rzędu ($1/w^2$).

**1. Pierwszy fragment (ze slajdu "The attractor – the transseries view z tej prezentacji na UJ "):**

$$\mathcal{A} = \frac{8C_\eta}{w} + \frac{16C_\eta C_{\tau_\pi}}{3w^2} + \dots$$

Ten wynik pochodzi z "czystej" **teorii Müllera-Israela-Stewarta (MIS)**. W tej teorii zakłada się najprostszą możliwą formę równania relaksacji dla naprężeń ścinających. Jest to najprostszy model, który jest przyczynowy i ma prawidłowy limit Naviera-Stokesa.

**2. Drugi fragment (z równania (7) w pracy prof Spalińskiego, "Initial State and Approach to Equilibrium"):**

$$ A(w) = \frac{8C_\eta}{w} + \frac{16C_\eta(C_{\tau_\pi} - C_{\lambda_1})}{3w^2} + O\left(\frac{1}{w^3}\right) $$

Ten wynik pochodzi z bardziej zaawansowanej **teorii BRSSS (Baier-Romatschke-Son-Starinets-Stephanov)**. Ta teoria jest uogólnieniem MIS. Autorzy BRSSS zdali sobie sprawę, że na poziomie drugiego rzędu w rozwinięciu gradientowym, oprócz członów obecnych w MIS, symetrie (kowariancja Lorentza i konforemność) dopuszczają istnienie dodatkowych członów. Najważniejszy z nich to człon sprzęgający tensor naprężeń ścinających sam ze sobą, którego siłę określa nowy współczynnik transportu, $\lambda_1$ (lub w notacji bezwymiarowej $C_{\lambda_1}$).

### Dlaczego Istnieje Wiele Teorii?

Hydrodynamika jest **teorią efektywną**, a nie fundamentalną. Oznacza to, że nie jest ona unikalna. Możemy skonstruować wiele różnych zestawów równań hydrodynamicznych, które:
a) Respektują fundamentalne prawa zachowania ($\nabla_\mu T^{\mu\nu}=0$).
b) Są przyczynowe (nie pozwalają na propagację sygnału z prędkością nadświetlną).
c) W granicy małych gradientów (blisko równowagi) zgadzają się z hydrodynamiką Naviera-Stokesa (czyli mają poprawny człon pierwszego rzędu, $\sim \eta/s$).

Teorie MIS i BRSSS są dwoma różnymi przykładami takich konstrukcji. Różnią się one tym, jak opisują **zachowanie systemu z dala od równowagi** (czyli poprawkami wyższych rzędów).

*   **Teoria MIS:** Jest najprostsza. Można ją postrzegać jako "minimalny model" przyczynowej hydrodynamiki.
*   **Teoria BRSSS:** Jest bardziej ogólna. Zawiera wszystkie możliwe człony do drugiego rzędu dozwolone przez symetrie. Posiada dodatkowe "pokrętła" (współczynniki transportu jak $C_{\lambda_1}$), które można dostroić tak, aby lepiej pasowały do wyników z fundamentalnej teorii mikroskopowej (np. z teorii kinetycznej lub holografii).

### Łopatologiczne Podsumowanie

Wyobraź sobie, że chcesz opisać ruch samochodu.

*   **Zerowy rząd (płyn idealny):** Mówisz, że samochód ma stałą prędkość. To bardzo proste, ale często nieprawdziwe.
*   **Pierwszy rząd (Navier-Stokes / MIS / BRSSS):** Uwzględniasz tarcie (opór powietrza). Mówisz: "siła oporu jest proporcjonalna do prędkości". To jest człon $\sim 1/w$. Wszystkie sensowne teorie (MIS, BRSSS) muszą się tu zgadzać, bo to jest najważniejszy efekt lepki. Dlatego współczynnik $a_1 = 8C_\eta$ jest **uniwersalny** dla wszystkich tych teorii.
*   **Drugi rząd (różnice między teoriami):** Teraz chcesz być bardziej dokładny.
    *   **Model MIS** mówi: "Wystarczy, uwzględniam tylko podstawowe tarcie". Jego rozwinięcie drugiego rzędu będzie miało jakąś prostą formę, np. $\frac{16C_\eta C_{\tau_\pi}}{3w^2}$.
    *   **Model BRSSS** mówi: "Chwileczkę, przy dużych prędkościach pojawiają się też turbulencje i inne efekty! Muszę dodać nowy człon opisujący, jak opór zależy od kwadratu prędkości". Ten nowy efekt jest opisany przez współczynnik $C_{\lambda_1}$. Dlatego w rozwinięciu BRSSS pojawia się dodatkowy kawałek: $- \frac{16C_\eta C_{\lambda_1}}{3w^2}$.

**Oba modele są "poprawne" jako teorie efektywne.** Który jest "lepszy"? Ten, który lepiej zgadza się z prawdziwym, mikroskopowym opisem ruchu samochodu (czyli z wynikami z tunelu aerodynamicznego). W fizyce QGP, rolę tunelu aerodynamicznego pełnią obliczenia w teorii kinetycznej lub holografii. Okazuje się, że te mikroskopowe teorie "preferują" opis BRSSS, ponieważ one również generują człon związany z $C_{\lambda_1}$.

Dlatego w nowszych pracach  używa się teorii BRSSS, a w starszych lub bardziej uproszczonych (jak na slajdzie, który może być bardziej pedagogiczny) używa się prostszej teorii MIS. **Różne rozwiązania odzwierciedlają po prostu użycie różnych modeli hydrodynamicznych do opisu tej samej fizyki.**
