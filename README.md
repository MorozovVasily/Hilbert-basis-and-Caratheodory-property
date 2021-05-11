# Hilbert-basis-and-Caratheodory-property
Код для выпускной квалификационной работы по теме "Базисы Гильберта и свойство Каратеодори"

Для запуска необходимо сначала установить библиотеку [Libnormalz](https://github.com/Normaliz/Normaliz) и ее зависимости.

Для компиляции я использовал следующую комманду:
```bash
g++-10 main.cpp -std=c++20 -O3 -I./pathToNormalizHeaders/ -L/pathToNormalizLibFiles -lnormaliz -lgmp -lgmpxx -leantic -lflint -larb -fopenmp -lpthread -Wall -Wno-sign-compare
```