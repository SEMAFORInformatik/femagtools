
    
class Problem(object):
    def __init__(self, n, ni, nf):
        self.dimension = n
        self.i_dimension = ni
        self.f_dim = nf
        self.c_dim = 0
        self.ic_dim = 0
        self.c_tol = []
        self.set_bounds([0.0]*n, [1.0]*n)

    def set_bounds(self, lb, ub):
        self.lower = lb
        self.upper = ub
        
    def test_constraint(self, c, i):
        if i < self.c_dim - self.ic_dim:
            return abs(c[i]) <= self.c_tol[i]
        return c[i] <= self.c_tol[i]

    def feasibility_c(self, c):
        for i, k in enumerate(c):
            if not self.test_constraint(c, i):
                return False
        return True
    
    def compare_constraints(self, c1, c2):
        count1 = count2 = 0
        norm1 = norm2 = 0
        # equality constraints
        for i in range(self.c_dim-self.ic_dim):
            if self.test_contraint(c1, i):
                count1 += 1
            if self.test_contraint(c2, i):
                count2 += 1
            norm1 += abs(c1[i]) * abs(c1[i])
            norm2 += abs(c2[i]) * abs(c2[i])
        # in-equality constraints
        for i in range(self.c_dim-self.ic_dim, self.c_dim):
            if self.test_contraint(c1, i):
                count1 += 1
            else:
                norm1 += c1[i] * c1[i]
            if self.test_contraint(c2, i):
                count2 += 1
            else:
                norm2 += c2[i] * c2[i]
        if count1 > count2:
            return True
        if count1 < count2:
            return False
        return norm1 < norm2
        
    def compare_fitness(self, v_f1, v_f2):
        count1 = count2 = 0
        for f1, f2 in zip(v_f1, v_f2):
            if f1 < f2:
                count1 += 1
            elif f1 == f2:
                count2 += 1
        return count1 + count2 == len(v_f1) and count1 > 0

    def compare_fc_impl(self, f1, c1, f2, c2):
        test1 = self.feasibility_c(c1)
        test2 = self.feasibility_c(c2)
        if test1 and not test2:
            return True
        if not test1 and test2:
            return False
        if test1:
            return self.compare_fitness(f1, f2)
        return self.compare_contraints(c1, c2)
    
    def compare_fc(self, f1, c1, f2, c2):
        if self.c_dim > 0:
            return self.compare_fc_impl(f1, c1, f2, c2)
        return self.compare_fitness(f1, f2)
        
