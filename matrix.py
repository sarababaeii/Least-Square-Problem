import sys


class Matrix:
    def __init__(self, rows, columns):
        self.rows = rows
        self.columns = columns
        self.entries = [[0 for x in range(columns)] for y in range(rows)]

    @staticmethod
    def zero(rows, columns):
        matrix = Matrix(rows, columns)
        return matrix

    def input(self):
        for i in range(self.rows):
            for j in range(self.columns):
                s = i.__str__() + ", " + j.__str__() + ": "
                self.entries[i][j] = float(input(s))

    def output(self):
        for i in range(self.rows):
            for j in range(self.columns):
                # s = i.__str__() + ", " + j.__str__() + ": "
                # print(s, end=": ")
                print(self.entries[i][j], end=" ")
            print("")

    def clone(self):
        matrix = Matrix(self.rows, self.columns)
        for i in range(self.rows):
            for j in range(self.columns):
                matrix.entries[i][j] = self.entries[i][j]
        return matrix

    def equals_to(self, matrix):
        for i in range(self.rows):
            for j in range(self.columns):
                if self.entries[i][j] != matrix.entries[i][j]:
                    return False
        return True

    def transpose(self):
        matrix = Matrix(self.columns, self.rows)
        for i in range(self.rows):
            for j in range(self.columns):
                matrix.entries[j][i] = self.entries[i][j]
        return matrix

    def multiply(self, matrix):
        ans = Matrix(self.rows, matrix.columns)
        for i in range(ans.rows):
            for j in range(ans.columns):
                for k in range(self.columns):
                    t = self.entries[i][k] * matrix.entries[k][j]
                    ans.entries[i][j] += t
        return ans

    def get_sub_matrix(self, from_row):
        matrix = SquareMatrix(self.rows - from_row)
        for i in range(matrix.rows):
            for j in range(matrix.columns):
                matrix.entries[i][j] = self.entries[i + from_row][j + from_row]
        return matrix

    def get_column_vector(self, column):
        v = Vector(self.rows)
        for i in range(self.rows):
            v.entries[i][0] = self.entries[i][column]
        return v

    def get_column_index_with_maximum_norm(self):
        max_norm = 0
        index = -1
        for i in range(self.columns):
            v = self.get_column_vector(i)
            norm = v.norm_2()
            if norm > max_norm:
                max_norm = norm
                index = i
        return index

    def qr_decomposition(self):
        square_matrix = SquareMatrix.squarization_by_insert(self)
        q, r = square_matrix.qr_decomposition_with_pivoting()
        r = r.desquarization(self.rows, self.columns)
        return [q, r]

    def augmented(self, b):
        matrix = Matrix(self.rows, self.columns + 1)
        for i in range(matrix.rows):
            for j in range(matrix.columns):
                if j < self.columns:
                    matrix.entries[i][j] = self.entries[i][j]
                else:
                    matrix.entries[i][j] = b.entries[i][0]
        return matrix

    def gaussian_elimination_method(self):
        a = self.__lu_decomposition()
        x = a.__backward_substitution()
        return x

    def __lu_decomposition(self):
        n = self.rows
        for i in range(n):
            if self.entries[i][i] == 0.0:
                sys.exit('Singular Matrix')
            for j in range(i + 1, n):
                ratio = self.entries[j][i] / self.entries[i][i]
                for k in range(n + 1):
                    self.entries[j][k] = self.entries[j][k] - ratio * self.entries[i][k]
        return self

    def __backward_substitution(self):
        n = self.rows
        x = Vector(n)
        x.entries[n - 1][0] = self.entries[n - 1][n] / self.entries[n - 1][n - 1]
        for i in range(n - 2, -1, -1):
            x.entries[i][0] = self.entries[i][n]
            for j in range(i + 1, n):
                x.entries[i][0] = x.entries[i][0] - self.entries[i][j] * x.entries[j][0]
            x.entries[i][0] = x.entries[i][0] / self.entries[i][i]
        return x

from vector import Vector
from square_matrix import SquareMatrix