from matrix import Matrix
from math import sqrt


class SquareMatrix(Matrix):
    def __init__(self, dimension):
        super(SquareMatrix, self).__init__(dimension, dimension)
        self.dimension = dimension

    @staticmethod
    def cast_to_square_matrix(matrix):
        square_matrix = SquareMatrix(matrix.rows)
        for i in range(square_matrix.dimension):
            for j in range(square_matrix.dimension):
                square_matrix.entries[i][j] = matrix.entries[i][j]
        return square_matrix

    @staticmethod
    def identity(dimension):
        matrix = SquareMatrix(dimension)
        for i in range(dimension):
            matrix.entries[i][i] = 1
        return matrix

    @staticmethod
    def squarization_by_insert(matrix):
        is_vertical = True
        square_matrix = SquareMatrix(matrix.rows)
        if matrix.rows < matrix.columns:
            square_matrix = SquareMatrix(matrix.columns)
            is_vertical = False
        for i in range(square_matrix.dimension):
            for j in range(square_matrix.dimension):
                if (is_vertical and j >= matrix.columns) or (not is_vertical and i >= matrix.rows):
                    square_matrix.entries[i][j] = 0
                else:
                    square_matrix.entries[i][j] = matrix.entries[i][j]
        return square_matrix

    @staticmethod
    def squarization_by_eliminate(matrix):
        square_matrix = SquareMatrix(matrix.columns)
        if matrix.rows < matrix.columns:
            square_matrix = SquareMatrix(matrix.rows)
        for i in range(square_matrix.dimension):
            for j in range(square_matrix.dimension):
                square_matrix.entries[i][j] = matrix.entries[i][j]
        return square_matrix

    def desquarization(self, row, column):
        matrix = Matrix(row, column)
        for i in range(row):
            for j in range(column):
                matrix.entries[i][j] = self.entries[i][j]
        return matrix

    def clone(self):
        matrix = SquareMatrix(self.dimension)
        for i in range(self.dimension):
            for j in range(self.dimension):
                matrix.entries[i][j] = self.entries[i][j]
        return matrix

    def qr_decomposition(self):
        q = SquareMatrix.identity(self.dimension)
        r = self.clone()
        for j in range(r.dimension):
            new_q, r = r.__make_column_zero(j)
            q = q.multiply(new_q)
            q = SquareMatrix.cast_to_square_matrix(q)
        return [q, r]

    def qr_decomposition_with_pivoting(self):
        q = SquareMatrix.identity(self.rows)
        r = self.clone()
        for j in range(self.dimension):
            if r.__can_continue_decomposition(j):
                r = r.__exchange_maximum_norm_column_with_first_column(j)  # permutation matrix?
                new_q, r = r.__make_column_zero(j)
                q = q.multiply(new_q)
                q = SquareMatrix.cast_to_square_matrix(q)
        return [q, r]

    def __can_continue_decomposition(self, column):
        remainder = self.get_sub_matrix(column)
        return not remainder.equals_to(Matrix.zero(remainder.rows, remainder.columns))

    def __make_column_zero(self, column):
        q = SquareMatrix.identity(self.dimension)
        r = self.clone()
        for i in range(column + 1, r.dimension):
            new_q, r = r.__make_entry_zero(i, column)
            q = q.multiply(new_q)
            q = SquareMatrix.cast_to_square_matrix(q)
        return [q, r]

    def __make_entry_zero(self, row, column):
        q = SquareMatrix.identity(self.dimension)
        r = self.clone()
        if r.entries[row][column] != 0:
            new_q = r.__rotator(row, column)
            q = q.multiply(new_q)
            q = SquareMatrix.cast_to_square_matrix(q)
            new_q_t = new_q.transpose()
            new_q_t = SquareMatrix.cast_to_square_matrix(new_q_t)
            r = new_q_t.multiply(r)
            r = SquareMatrix.cast_to_square_matrix(r)
        return [q, r]

    def __rotator(self, row, column):
        a = self.entries[column][column]
        b = self.entries[row][column]
        length = sqrt((a * a) + (b * b))
        sin = b / length
        cos = a / length

        q = SquareMatrix.identity(self.dimension)
        q.entries[column][column] = cos
        q.entries[column][row] = -sin
        q.entries[row][column] = sin
        q.entries[row][row] = cos
        return q

    def __exchange_maximum_norm_column_with_first_column(self, index):
        matrix = self.clone()
        sub_matrix = matrix.get_sub_matrix(index)
        maximum_column_index = sub_matrix.get_column_index_with_maximum_norm()
        if maximum_column_index != 0:
            matrix = matrix.__exchange_columns(index, index + maximum_column_index)
        return matrix

    def __exchange_columns(self, first, second):
        matrix = SquareMatrix(self.dimension)
        for i in range(self.dimension):
            for j in range(self.dimension):
                if j == first:
                    matrix.entries[i][j] = self.entries[i][second]
                elif j == second:
                    matrix.entries[i][j] = self.entries[i][first]
                else:
                    matrix.entries[i][j] = self.entries[i][j]
        return matrix
