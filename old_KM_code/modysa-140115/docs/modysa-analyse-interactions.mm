= MoDySA =

Opis pakietu do analizy trajektorii atomowych uzyskanych na drodze symulacji dynamiki molekularnej.

== modysa.analyse.interactions ==

 * Opis konwencji w nazewnictwie grup atomów definiowanych na potrzeby analizy oddziaływań w obrębie interfazy woda/dwuwarstwa lipidowa.
[[Anchor(ndx_file)]]

 W celu przeprowadzenia analizy wiązań wodorowych oraz pomostów wodnych należy utworzyć plik indeksowy zawierający definicję następujących grup:

   donor_atoms:: indeksy donorów wiązania wodorowego. Indeksy atomów z tej samej cząsteczki powinny być podane w tym samym wierszu. Zauważmy, że taka definicja formatu grupy umożliwia umieszczenie w `donor_atoms` indeksów atomów z różnych (z chemicznego punktu widzenia) cząsteczek. Istnienie tej i kolejnej grupy '''ma umożliwić''' rozróżnianie w obliczeniach oddziaływań wewnątrz- i międzycząsteczkowych.

   np. zapis:

   {{{[ donor_atoms ]
1 3 8  # cząsteczka POPE nr 1
3224 3229  # cząsteczka POPG nr 8
}}}

   oznacza, że mamy zdefinowane pięć atomów donorowych w dwóch cząsteczkach. ``Walidacja danych``: indeksy nie mogą się powtarzać, największy indeks dla danej cząsteczki nie może być większy niż najmniejszy indeks dla innej cząsteczki. Istnieje możliwość umieszczania komentarzy na końcu każdego wiersza po znaku #. Inną zaletą takiego formatu jest możliwość definiowania fragmentów cząsteczek lub reszt (wystarczy odpowiednie atomy umieścić w jednej linii) dzięki czemu można później analizować np. oddziaływania w obrębie reszty, wzgl. między resztami.

   acceptor_atoms:: indeksy akceptorów wiązania wodorowego. Indeksy atomów z tej samej cząsteczki powinny być podane w tym samym wierszu.

   H_count:: liczba atomów wodoru związanych z ''i''-tym atomem donora. ''Walidacja danych'': (1) liczba atomów w `donor_atoms` musi być równa liczbie atomów w `H_count`, (2) suma liczb w tej grupie powinna być równa liczbie atomów w `H_atoms`, (3) liczba atomów wodoru musi wynosić co najmniej 1. np. zapis:

   {{{[ H_count ]
1 3 2
1 1
}}}

   oznacza, że atomy pierwszej ze zdefinowanych cząsteczek (POPE nr 1) są związane kolejno z: 1, 3 i 2 atomami wodoru

   H_atoms:: indeksy atomów wodoru związanych z kolejnymi atomami donorów. Indeksy atomów z tej samej cząsteczki powinny być podane w tym samym wierszu. ''Walidacja danych'': (1) suma atomów wodoru zdefinowana dla poszczególnych cząsteczek powinna być równa sumie liczb z odpowiedniego wiersza grupy `H_count`, (2) liczba cząsteczek zdefinowanych w grupach `donor_atoms`, `H_count` i `H_atoms` powinna być taka sama. Przykład:

   {{{[ H_atoms ]
2   4 5 6   9 10
3225   3230
}}}

   molecules:: indeksy do grup `donor_count` oraz `acceptor_count` wskazujące na atomy (donorów lub akceptorów) należące do tej samej cząsteczki. Jeśli dana cząsteczka zawiera wyłącznie atomy donorowe lub wyłącznie atomy akceptorowe, wtedy na pierwszej bądź drugiej pozycji w parze liczb przypisanych poszczególnym cząsteczkom podana jest wartość `-1`. Zauważmy, że nawet jeśli cząsteczka zawiera atomy donorowe i akceptorowe, możemy tu wymusić wyznaczanie oddziaływań tylko określonego typu: np. obejmujące tylko takie wiązania wodorowe, w których donorem jest cząsteczka 1 a akceptor jest w innych wskazanych (patrz grupa `pairs`) cząsteczkach. Konsekwentnie indeksy cząsteczek w `donor_count` i `acceptor_count` podajemy dla każdej definowanej cząsteczki w jednym wierszu. ''Walidacja danych'': (1) liczba indeksów powinna być parzysta ponieważ dla każdej cząsteczki mamy zawsze parę liczb, indeksy na pierwszych pozycjach dwuelementowych list przypisanych kolejnym cząsteczkom powinny być z zakresu od 1 do ''n'' i od 1 do ''k'', gdzie ''n'' jest liczbą elementów grupy `donor_count` a ''k'' - liczbą elementów grupy `acceptor_count`; (2) Nie jest dopusczalne wystąpienie wiersza `-1 -1` ponieważ oznaczałoby to, że definiujemy cząsteczkę, która nie ma ani donorów ani akceptorów. Zgodnie z opisanym formatem możemy dla fizycznie jednej cząsteczki zdefinować max. trzy warianty, np. `1 1`, `1 -1` oraz `-1 1`. Na obecną chwilę nie widzę zastosowań takiej definicji jednak warto pamiętać,że taka możliwość istnieje. Przykład:

   {{{[ molecules ]
1 1
2 -1  # wariant POPG:exclusive_donor
-1 2  # wariant POPG:exclusive_acceptor
2 2   # wariant POPG:donor_and_acceptor
}}}

   pairs:: para indeksów cząsteczek, dla których chcemy wyznaczyć wiązania wodorowe lub pomosty wodne. Indeksy odnoszą się do pozycji w grupie `molecules` i wskazują zawsze numer porządkowy cząsteczki (od 1 do ''n''). Zauważmy, że można w tym miejscu kontrolować sposób wyznaczania oddziaływań, tj. między- czy wewnątrzcząsteczkowe, np. zapis: 

   {{{[ pairs ]
1 1   
1 2 
1 3   
3 3
}}}

   oznacza, że interesują nas oddziaływania wewnątrz-cząsteczkowe dla cząsteczki nr 1 i 3, a także oddziaływania międzycząsteczkowe dla cząsteczek 1 i 2 oraz 1 i 3. ''Walidacja danych'': indeksy muszą wskazywać na istniejące cząsteczki, tj. np. jeśli największy indeks cząsteczki w `molecules` wynosi `n`, to indeks w `pairs` musi być niewiększy niż `n`.

-----

 * w module `modysa.analyse.interactions` lub gdzie indziej (może: `modysa.analyse.aux`) powinny być zdefiniowane funkcje pomocnicze, umożliwiające tworzenie odpowiednich plików indeksowych. Ponieważ generowanie plików indeksowych może być w ogólnym przypadku dosyć skomplikowane być może warto obsługę tego zadania przenieść do wydzielonego modułu (np. `modysa.io.ndx`)
