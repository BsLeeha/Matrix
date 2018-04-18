# Matrix 介绍
**主要功能**: </br>
基本矩阵运算的实现，基本数值分析方法的实现

**所用方法**: </br>
C++ template

**参考资料**: </br>
David C.Lay *线性代数及其应用* </br>
Timothy Sauer *数值分析* </br>
其它(主要来自 google) 

# 各头文件主要功能介绍

**Polynomial.hpp**: </br>
Poly 类的定义

**Matrix.hpp**: </br>
Matrix 类的定义及四则运算等实用函数

**Basic.hpp**: </br>
矩阵的一些分割；</br>
特殊矩阵的生成：单位矩阵，随机矩阵； </br>
矩阵的基本计算：求幂，求逆，置换，转置，求秩，行列式，2范数，最小二乘； </br>
等

**Elimination.hpp**：</br>
矩阵消元：求主元，上三角，下三角，行最简阶梯形；

**Factorization.hpp**: </br>
矩阵分解：PLU 分解，QR 分解的几种实现(Gram-Scmidt);

**EigenValueEstimate**: </br>
特征值和特征向量的求取：幂法，反幂法；

**SystemSolving**: </br>
方程组求解：直接法(两种高斯消元)，迭代法(Jacobi 迭代，GaussSeidel 迭代);


