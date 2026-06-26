
import math


def evaluate(x_in, params=[]):
    _, sample = x_in
    u_1 = sample['u_1']
    u_2 = sample['u_2']
    x = math.exp(sample['x'])
    output = 2. + 0.4 * u_1 - 0.2 * u_2 + 0.3 * x + 0.01 * x * u_1
    return (output,)
