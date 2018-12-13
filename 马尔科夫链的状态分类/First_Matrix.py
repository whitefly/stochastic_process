"""
本代码是为了解决随机过程的一道习题(例4.2.5): 根据转移矩阵求n阶首达矩阵 的暴力计算型解法
基于符号运算库sympy
核心思想:n阶的首达矩阵需要:先求出1~n阶转移矩阵,1~n-1阶首达矩阵,然后根据公式求解
"""
from sympy import symbols, eye, Matrix, simplify


class Solver:
    def __init__(self, P):
        self._max = 100
        self.P = P
        self._PM_size = 1  # 已经计算出来的转移矩阵的阶数
        self._fm_size = 1  # 已经计算出来的首达矩阵的阶数
        self._P_M = [0, P] + [0] * (self._max - 2)  # 转移概率矩阵,P_M[1]就是1步转移矩阵
        self._f_M = [0, P] + [0] * (self._max - 2)  # 首达矩阵,f_M[1]就是1步首达概率矩阵

    def _generate_f_m1(self, n):
        # 生成n阶首达概率矩阵
        m = self.P.shape[0]
        temp = eye(m, m)
        # 根据公式计算出每个元素值
        for i in range(m):
            for j in range(m):
                pij_n = self._P_M[n][i, j]
                total = 0
                for L in range(1, n):
                    total += self._f_M[L][i, j] * self._P_M[n - L][j, j]
                temp[i, j] = simplify(pij_n - total)
        return temp

    def _generate_p_m(self, n):
        # 生成n阶转移矩阵
        if self._PM_size >= n:
            raise ValueError('{}阶转移矩阵已经存在'.format(n))
        for i in range(self._PM_size + 1, n + 1):
            self._P_M[i] = self._P_M[i - 1] * self.P

        self._PM_size = n

    def _generate_f_m(self, n):
        # 为生成n阶的首达矩阵创造条件:求n阶首达矩阵,需要先求出1~n阶转移矩阵,1~n-1阶首达矩阵
        if self._PM_size < n:
            self._generate_p_m(n)

        if self._fm_size < n:
            for i in range(self._fm_size + 1, n + 1):
                self._f_M[i] = self._generate_f_m1(i)

        self._fm_size = n

    def get_fm(self, n):
        if self._fm_size < n:
            self._generate_f_m(n)
        return self._f_M[n]


if __name__ == '__main__':
    p1, q1, q2, p2, p3, q3 = symbols('p1,q1,q2,p2,p3,q3')
    P = Matrix([[0, p1, q1], [q2, 0, p2], [p3, q3, 0]])
    s = Solver(P)
    print(s.get_fm(4))
