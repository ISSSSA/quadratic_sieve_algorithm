from itertools import chain
from math import sqrt, exp, log

from src import STonelli, is_probable_prime


def gcd(a, b):  # Алгоритм Евклида
    if b == 0:
        return a
    elif a >= b:
        return gcd(b, a % b)
    else:
        return gcd(b, a)


def isqrt(n):  # Метод Ньютона
    x = n
    y = (x + 1) // 2
    while y < x:
        x = y
        y = (x + n // x) // 2
    return x


def mprint(M):  # принт матрицы
    for row in M:
        print(row)


def prime_gen(n):  # генерация числе решетом Эратосфена
    if n < 2:
        return []

    nums = []
    isPrime = []

    for i in range(0, n + 1):
        nums.append(i)
        isPrime.append(True)

    isPrime[0] = False
    isPrime[1] = False

    for j in range(2, int(n / 2)):
        if isPrime[j] == True:
            for i in range(2 * j, n + 1, j):
                isPrime[i] = False

    primes = []
    for i in range(0, n + 1):
        if isPrime[i] == True:
            primes.append(nums[i])

    return primes


def quad_residue(a, n):
    # Квадратичный остаток
    l = 1
    q = (n - 1) // 2
    x = q ** l
    if x == 0:
        return 1

    a = a % n
    z = 1
    while x != 0:
        if x % 2 == 0:
            a = (a ** 2) % n
            x //= 2
        else:
            x -= 1
            z = (z * a) % n

    return z


def size_bound(N):  # ищет оптимальную базу для факторизации и размер

    F = pow(exp(sqrt(log(N) * log(log(N)))), sqrt(2) / 4)
    I = F ** 3
    # print(F,I)
    return int(F), int(I)


def find_base(N, B):
    # генерирует B

    factor_base = []
    primes = prime_gen(B)
    # print(primes)

    for p in primes:
        if quad_residue(N, p) == 1:
            factor_base.append(p)
    return factor_base


def find_smooth(factor_base, N, I):
    # используя решето ищет гладки-B числа

    def sieve_prep(N, sieve_int):
        sieve_seq = [x ** 2 - N for x in range(root, root + sieve_int)]
        return sieve_seq

    sieve_seq = sieve_prep(N, I)
    sieve_list = sieve_seq.copy()
    if factor_base[0] == 2:
        i = 0
        while sieve_list[i] % 2 != 0:
            i += 1
        for j in range(i, len(sieve_list), 2):
            while sieve_list[j] % 2 == 0:
                sieve_list[j] //= 2
    # print("")
    for p in factor_base[1:]:
        residues = STonelli(N, p)

        for r in residues:
            for i in range((r - root) % p, len(sieve_list), p):
                while sieve_list[i] % p == 0:
                    sieve_list[i] //= p
    xlist = []
    smooth_nums = []
    indices = []

    for i in range(len(sieve_list)):
        if len(smooth_nums) >= len(factor_base) + T:
            break
        if sieve_list[i] == 1:
            smooth_nums.append(sieve_seq[i])
            xlist.append(i + root)
            indices.append(i)

    return (smooth_nums, xlist, indices)


def build_matrix(smooth_nums, factor_base):
    # генерирует экспоненциальные векторы по модулю 2 из ранее полученных гладких чисел, затем строит матрицу

    def factor(n, factor_base):
        factors = []
        if n < 0:
            factors.append(-1)
        for p in factor_base:
            if p == -1:
                pass
            else:
                while n % p == 0:
                    factors.append(p)
                    n //= p
        return factors

    M = []
    factor_base.insert(0, -1)
    for n in smooth_nums:
        exp_vector = [0] * (len(factor_base))
        n_factors = factor(n, factor_base)
        for i in range(len(factor_base)):
            if factor_base[i] in n_factors:
                exp_vector[i] = (exp_vector[i] + n_factors.count(factor_base[i])) % 2

        if 1 not in exp_vector:  # поиск квадратов
            return True, n
        else:
            pass

        M.append(exp_vector)
    return (False, transpose(M))


def transpose(matrix):
    # транспонирование матрицы
    new_matrix = []
    for i in range(len(matrix[0])):
        new_row = []
        for row in matrix:
            new_row.append(row[i])
        new_matrix.append(new_row)
    return (new_matrix)


def gauss_elim(M):
    # сокращенная форма исключения по Гауссу, находит rref и считывает нулевое пространство

    marks = [False] * len(M[0])

    for i in range(len(M)):
        row = M[i]

        for num in row:
            if num == 1:
                j = row.index(num)
                marks[j] = True

                for k in chain(range(0, i), range(i + 1, len(M))):
                    if M[k][j] == 1:
                        for i in range(len(M[k])):
                            M[k][i] = (M[k][i] + row[i]) % 2
                break

    M = transpose(M)

    sol_rows = []
    for i in range(len(marks)):
        if marks[i] == False:
            free_row = [M[i], i]
            sol_rows.append(free_row)

    if not sol_rows:
        return ("Не найдено достаточно решений. Нужно больше гладких-чисел")
    print("Найдены {} возможные решения".format(len(sol_rows)))
    return sol_rows, marks, M


def solve_row(sol_rows, M, marks, K=0):
    solution_vec, indices = [], []
    free_row = sol_rows[K][0]
    for i in range(len(free_row)):
        if free_row[i] == 1:
            indices.append(i)
    for r in range(len(M)):
        for i in indices:
            if M[r][i] == 1 and marks[r]:
                solution_vec.append(r)
                break

    solution_vec.append(sol_rows[K][1])
    return (solution_vec)


def solve(solution_vec, smooth_nums, xlist, N):
    solution_nums = [smooth_nums[i] for i in solution_vec]
    x_nums = [xlist[i] for i in solution_vec]

    Asquare = 1
    for n in solution_nums:
        Asquare *= n

    b = 1
    for n in x_nums:
        b *= n

    a = isqrt(Asquare)

    factor = gcd(b - a, N)
    return factor


def QS(n, B, I):
    # однополиномиальная версия квадратичного сита с заданной границей гладкости B и интервалом сита I

    global N
    global root
    global T
    N, root, K, T = n, int(sqrt(n)), 0, 1

    if is_probable_prime(N):
        return "простое"

    if isinstance(sqrt(N), int):
        return isqrt(N)

    print("Попытка фактаризовать {}...".format(N))

    print("Генерация {}-гладкочисленной базы...".format(B))
    factor_base = find_base(N, B)

    global F
    F = len(factor_base)

    print("Поиск {} {}-гладких отношений...".format(F + T, B))
    smooth_nums, xlist, indices = find_smooth(factor_base, N, I)

    print("Найдены {} B-гладкие числа.".format(len(smooth_nums)))

    print(smooth_nums)

    if len(smooth_nums) < len(factor_base):
        return ("Не достаточно гладких-чисел. Увеличте интервал решета или размер границ.")

    print("Построение экспоненциальной матрицы...")
    is_square, t_matrix = build_matrix(smooth_nums, factor_base)
    # builds exponent matrix mod 2 from relations

    if is_square == True:
        x = smooth_nums.index(t_matrix)
        factor = gcd(xlist[x] + sqrt(t_matrix), N)
        print("Квадрат найден!")
        return factor, N / factor

    print("Выполняем исключения Гаусса...")
    sol_rows, marks, M = gauss_elim(t_matrix)
    solution_vec = solve_row(sol_rows, M, marks, 0)

    print("Решаем квадратичное сравнение...")

    factor = solve(solution_vec, smooth_nums, xlist, N)

    for K in range(1, len(sol_rows)):
        if (factor == 1 or factor == N):
            print("Пробуем другой вектор для решения...")
            solution_vec = solve_row(sol_rows, M, marks, K)
            factor = solve(solution_vec, smooth_nums, xlist, N)
        else:
            print("Найдены простые множители!")
            return factor, int(N / factor)

    return ("Не найдены нетриваильные простые множители")


if __name__ == "__main__":
    number_to_factorize = int(input("Введите число которое хотите факторизовать: "))
    B = int(input("Введите границу гладкости B: "))
    I = int(input("Введите интервал решета I: "))
    print(QS(number_to_factorize, B, I))
