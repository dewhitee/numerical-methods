![all labworks](https://media.discordapp.net/attachments/578480014740422676/783811352182652938/unknown.png)
#### 1 лаб 
  Метод исключения Гаусса с ведущим элементом;
  
  Метод Гаусса-Зейделя;
  
  Эесперементальное определение числа обусловленности матрицы;
#### 2 лаб   
  Аппроксимация по методу наименьших квадратов(МНК);
 
  Сплайн интерполяция;
  
  Графическое сравнение двух методов;
#### 3 лаб  
  Метод бисекции(деление отрезка еа пополам);
  
  Метод парабол;
#### 4 лаб
  Метод Эйлера;
  
  Метод Рунге-Кутта 4-го порядка;
# 1 ЛАБ
1) Гаусс с ведущим элементом 

![alt text](https://github.com/dewhitee/numerical-methods/blob/main/Images/Gauss%20i%20vedusij%20element.png)

2)Гаусс-Зейдель

![alt text](https://github.com/dewhitee/numerical-methods/blob/main/Images/Method%20Gaussa%20Zeidela.png)

3)

Показывает точность потенциальной обусловленности 
Число обусловленности матрицы - это индикатор точности решения системы

### Норма вектора - ||x||
### Евклидова норма : ||x||=sqrt(x1^2 + x2^2 + x3^2 + .....)
### Cond(A)- коэффициент усиления ошибки 
По хорошему cond(A) 1 до 10 - это идеально
от 1 - 1000 - терпимо
больше 1000 - всё плохо, я в агонии

![alt text](https://media.discordapp.net/attachments/578480014740422676/783643689959424000/unknown.png)

![alt text](https://github.com/dewhitee/numerical-methods/blob/main/Images/%D0%A7%D0%B8%D1%81%D0%BB%D0%BE%20%D0%BE%D0%B1%D1%83%D1%81%D0%BB%D0%BE%D0%B2%D0%BB%D0%B5%D0%BD%D0%BD%D0%BE%D1%81%D1%82%D0%B8%20%D0%BC%D0%B0%D1%82%D1%80%D0%B8%D1%86%D1%8B.png?raw=true)

#### det(A) - детерминант (определитель) (прос)
det(A)=0    если    cond(A) = ∞ (бессконечность)

![alt text](https://wikimedia.org/api/rest_v1/media/math/render/svg/5b2e40d390e1d26039aabee44c7d1d86c8755232)
![alt text](https://wikimedia.org/api/rest_v1/media/math/render/svg/a891ca1b518ba39ff21a458c74f9cc74bcefb18c)

### Норма матрицы - ||A||=max||aj||
Матрицу можно рассматривать как вектор векторов, поэтому принято считать нормой матрицы - максимальную норму его столбцов.
![alt text](https://media.discordapp.net/attachments/578480014740422676/783348864052363284/unknown.png)

### Эксперементальное определение: 
![alt text](https://github.com/dewhitee/numerical-methods/blob/main/Images/Эксперементальный%20метод.png)


### Виды матриц:

![alt text](https://github.com/dewhitee/numerical-methods/blob/main/Images/Oboznacenie%20matrici%201.png)
![alt text](https://github.com/dewhitee/numerical-methods/blob/main/Images/Oboznacenie%20matrici%202.png)

# 2 лаб   
#### 1)Аппроксимация по методу наименьших квадратов(МНК);
 
  ![alt text](https://github.com/dewhitee/numerical-methods/blob/main/Images/%D0%90%D0%BF%D0%BF%D1%80%D0%BE%D0%BA%D1%81%D0%B8%D0%BC%D0%B0%D1%86%D0%B8%D1%8F.png)
  ![alt text](https://github.com/dewhitee/numerical-methods/blob/main/Images/%D0%9C%D0%B5%D1%82%D0%BE%D0%B4%20%D0%9C%D0%9D%D0%9A.png)
  ![alt text](https://github.com/dewhitee/numerical-methods/blob/main/Images/%D0%9C%D0%B5%D1%82%D0%BE%D0%B4%20%D0%9C%D0%9D%D0%9A%202.png)
  
#### 2)Сплайн интерполяция
  
 ![alt text](https://github.com/dewhitee/numerical-methods/blob/main/Images/%D0%9A%D1%83%D0%B1%D0%B8%D1%87%D0%B5%D1%81%D0%BA%D0%B8%D0%B9%20%D1%81%D0%BF%D0%BB%D0%B0%D0%B9%D0%BD.png)
 ![alt text](https://github.com/dewhitee/numerical-methods/blob/main/Images/%D0%9A%D1%83%D0%B1%D0%B8%D1%87%D0%B5%D1%81%D0%BA%D0%B8%D0%B9%20%D1%81%D0%BF%D0%BB%D0%B0%D0%B9%D0%BD%20%D0%B0%D0%BB%D0%B3%D0%BE%D1%80%D0%B8%D1%82%D0%BC%202.png)
 ![alt text](https://github.com/dewhitee/numerical-methods/blob/main/Images/%D0%9A%D1%83%D0%B1%D0%B8%D1%87%D0%B5%D1%81%D0%BA%D0%B8%D0%B9%20%D1%81%D0%BF%D0%BB%D0%B0%D0%B9%D0%BD%20%D0%B0%D0%BB%D0%B3%D0%BE%D1%80%D0%B8%D1%82%D0%BC%203.png)
 ![alt text](https://github.com/dewhitee/numerical-methods/blob/main/Images/%D0%9A%D1%83%D0%B1%D0%B8%D1%87%D0%B5%D1%81%D0%BA%D0%B8%D0%B9%20%D1%81%D0%BF%D0%BB%D0%B0%D0%B9%D0%BD%204.png)
 
 ##### Кубический сплайн алгоритм:
 
 ![alt text](https://github.com/dewhitee/numerical-methods/blob/main/Images/%D0%9A%D1%83%D0%B1%D0%B8%D1%87%D0%B5%D1%81%D0%BA%D0%B8%D0%B9%20%D1%81%D0%BF%D0%BB%D0%B0%D0%B9%D0%BD%20%D0%B0%D0%BB%D0%B3%D0%BE%D1%80%D0%B8%D1%82%D0%BC.png)
 ![alt text](https://github.com/dewhitee/numerical-methods/blob/main/Images/%D0%9A%D1%83%D0%B1%D0%B8%D1%87%D0%B5%D1%81%D0%BA%D0%B8%D0%B9%20%D1%81%D0%BF%D0%BB%D0%B0%D0%B9%D0%BD%20%D0%B0%D0%BB%D0%B3%D0%BE%D1%80%D0%B8%D1%82%D0%BC%202.png)
 ![alt text](https://github.com/dewhitee/numerical-methods/blob/main/Images/%D0%9A%D1%83%D0%B1%D0%B8%D1%87%D0%B5%D1%81%D0%BA%D0%B8%D0%B9%20%D1%81%D0%BF%D0%BB%D0%B0%D0%B9%D0%BD%20%D0%B0%D0%BB%D0%B3%D0%BE%D1%80%D0%B8%D1%82%D0%BC%203.png)
 ![alt text](https://github.com/dewhitee/numerical-methods/blob/main/Images/%D0%9A%D1%83%D0%B1%D0%B8%D1%87%D0%B5%D1%81%D0%BA%D0%B8%D0%B9%20%D1%81%D0%BF%D0%BB%D0%B0%D0%B9%D0%BD%20%D0%B0%D0%BB%D0%B3%D0%BE%D1%80%D0%B8%D1%82%D0%BC%204.png)
  
  ##### 3)Графическое сравнение двух методов;
