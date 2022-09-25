from matrix import Matrix
from math import sqrt


class Vector(Matrix):
    def __init__(self, dimension):
        super(Vector, self).__init__(dimension, 1)
        self.dimension = dimension

    @staticmethod
    def cast_to_vector(matrix):
        v = Vector(matrix.rows)
        for i in range(v.dimension):
            v.entries[i][0] = matrix.entries[i][0]
        return v

    @staticmethod
    def insert_to_vector(vector, dimension):
        v = Vector(dimension)
        for i in range(vector.dimension):
            v.entries[i][0] = vector.entries[i][0]
        return v

    def norm_2(self):
        ans = 0
        for i in range(self.dimension):
            t = self.entries[i][0]
            ans += (t * t)
        return sqrt(ans)

    def block(self, rows):
        v = Vector(rows)
        for i in range(rows):
            v.entries[i][0] = self.entries[i][0]
        return v
