from matrix import Matrix
from vector import Vector
from square_matrix import SquareMatrix


def get_inputs():
    n = int(input("n: "))
    m = int(input("m: "))
    print("Enter Matrix A:")
    a = Matrix(n, m)
    a.input()
    print("Enter Vector b:")
    b = Vector(n)
    b.input()
    return [a, b]


def solve_least_square_problem(a, b):
    q, r = a.qr_decomposition()
    c = q.transpose().multiply(b)
    r_hat = SquareMatrix.squarization_by_eliminate(r)
    c_hat = Vector.cast_to_vector(c).block(r_hat.dimension)
    x_hat = solve_linear_equation_system(r_hat, c_hat)
    x = Vector.insert_to_vector(x_hat, b.dimension)
    return x


def solve_linear_equation_system(a, b):
    augmented_matrix = a.augmented(b)
    x = augmented_matrix.gaussian_elimination_method()
    return x


if __name__ == '__main__':
    a, b = get_inputs()
    x = solve_least_square_problem(a, b)
    print("x is:")
    x.output()
