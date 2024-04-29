import time

def BCD_IQ(x, h, U, D):
    # Start timing
    start = time.time()

    # Convolution of up_x and h
    y_long_size = len(x) * U + len(h) - 1
    y_size = round((y_long_size / D) + 0.5)
    y = [0] * y_size
    
    # Outer for loop ensures only every D-th element is calculated
    for n in range(0, y_long_size, D):
        # SHOULD NOT LOOP THROUGH EVERY h
        # Loop through every U h values, starting from n % U
        for k in range(n % U, len(h), U):
            if n-k >= 0 and (n-k) < (len(x) * U):
                y[int(n / D)] += h[k] * x[int((n-k) / U)]
    
    # End timing
    end = time.time()
    elapsed = end - start
    #print("Time for resample2: " + str(elapsed) + " seconds")
    return y, elapsed

def BCD_IQ(x, h, state, d):
    x_size = len(x)
    h_size = len(h)
    state_size = len(state)
    
    # Clear output vector
    y = [0] * (x_size // d)
    y_size = len(y)
    
    # Fast convolution/downsampling in one
    for i in range(y_size):
        for j in range(h_size):
            index = i * d - j
            if 0 <= index < x_size:
                y[i] += h[j] * x[index]
            elif index < 0:
                y[i] += h[j] * state[state_size + index]

    # Clear state vector (comment for improved performance?)
    # state.clear()  # You may skip this line if you want to keep the state from previous calls.
    # state.extend([0.0] * (len(h) - 1))

    # Copy end of input into state vector for next block
    state_index = 0
    for b in range(x_size - state_size, x_size):
        state[state_index] = x[b]
        state_index += 1
    return y, state


def deinterleave(in_list):
    size = len(in_list) // 2
    a = [0] * size
    b = [0] * size
    for i in range(size):
        a[i] = in_list[i * 2]
        b[i] = in_list[i * 2 + 1]
    return a, b



def list_avg(list):
    return sum(list) / len(list)

# Example usage
x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
h = [0.5, 1, 0.5]
U = 1
D = 2
#print("x length = " + str(len(x)))
#print("h length = " + str(len(h)))

# Output the resulting signal
#print("Resulting signal after fast resampling:")
#print(y_fast)
#print("Resulting signal after slow resampling:")
#print(y_slow)

flag = True
if len(y_fast) == len(y_slow):
    for i in range(len(y_fast)):
        if y_fast[i] != y_slow[i]:
            flag = False
            print("y_fast[" + str(i) + "] != y_slow[" + str(i))
            print("     y_fast = " + str(y_fast[i]) + ", y_slow = " + str(y_slow[i]))
    if flag == True:
        print("y_fast is equal to y_slow")
else:
    print("y_fast is length " + str(len(y_fast)) + ", y_slow is length " + str(len(y_slow)))