function x_cross = SkewSymmetric(x)

x_cross = [0, x(3), -x(2);
           -x(3), 0, x(1);
           x(2), -x(1), 0];