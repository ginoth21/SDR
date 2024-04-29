import time

def resample1(x, h, U, D):
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

def resample2(x, h, U, D):
    # Start timing
    start = time.time()

    # Convolution of up_x and h
    y_long_size = len(x) * U + len(h) - 1
    y_size = round((y_long_size / D) + 0.5)
    y = [0] * y_size
    for n in range(y_long_size):
        # This if could be incorporated in a for loop in C++ (n+=D)
        if(n % D == 0):
            for k in range(len(h)):
                if n-k >= 0 and (n-k) < (len(x) * U) and (n-k) % U == 0:
                    y[int(n / D)] += h[k] * x[int((n-k) / U)]
    
    # End timing
    end = time.time()
    elapsed = end - start
    #print("Time for resample2: " + str(elapsed) + " seconds")
    return y, elapsed

def slow_resample(x, h, U, D):
    # Start timing
    start = time.time()
    
    # Print original x signal
    #print("Original x:")
    #print(x)
    #print("")

    # Upsampling x signal
    up_x_size = len(x) * U
    up_x = [0] * up_x_size
    for i in range(len(x)):
        up_x[i * U] = x[i]
        # up_x[i] is only non-zero when i is a multiple of U
    #print("Upsampled x:")
    #print(up_x)
    #print("")

    # Convolution of up_x and h
    y_long_size = len(up_x) + len(h) - 1
    y_long = [0] * y_long_size
    for n in range(len(y_long)):
        for k in range(len(h)):
            if n-k >= 0 and n-k < len(up_x):
                # It's 0 whenever n-k is not a multiple of U
                y_long[n] += h[k] * up_x[n-k]
    #print("")
    #print("Long y:")
    #print(y_long)
    #print("")

    # Downsample y
    y_size = round((y_long_size / D) + 0.5)
    y = [0] * y_size
    for index in range(len(y)):
        y[index] = y_long[index * D]
        # Throws out y values where index is not a multiple of D

    # End timing
    end = time.time()
    elapsed = end - start
    #print("Time for slow_resample: " + str(elapsed) + " seconds")
    
    return y, elapsed

def list_avg(list):
    return sum(list) / len(list)

# Example usage
x = [1, 2, 3, 4, 5, 4, 6, 8, 2, 1, 1, 9, 10, 3, 7, 8, 1, 2, 3, 4, 5, 4, 6, 8, 2, 1, 1, 9, 10, 3, 7, 8, 1, 2, 3, 4, 5, 4, 6, 8, 2, 1, 1, 9, 10, 3, 7, 8, 1, 2, 3, 4, 5, 4, 6, 8, 2, 1, 1, 9, 10, 3, 7, 8, 1, 2, 3, 4, 5, 4, 6, 8, 2, 1, 1, 9, 10, 3, 7, 8, 1, 2, 3, 4, 5, 4, 6, 8, 2, 1, 1, 9, 10, 3, 7, 8, 1, 2, 3, 4, 5, 4, 6, 8, 2, 1, 1, 9, 10, 3, 7, 8, 1, 2, 3, 4, 5, 4, 6, 8, 2, 1, 1, 9, 10, 3, 7, 8, 1, 2, 3, 4, 5, 4, 6, 8, 2, 1, 1, 9, 10, 3, 7, 8, 1, 2, 3, 4, 5, 4, 6, 8, 2, 1, 1, 9, 10, 3, 7, 8, 1, 2, 3, 4, 5, 4, 6, 8, 2, 1, 1, 9, 10, 3, 7, 8, 1, 2, 3, 4, 5, 4, 6, 8, 2, 1, 1, 9, 10, 3, 7, 8]
h = [0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0.5, 0.5, 1, 0.5]
U = 441
D = 1600
#print("x length = " + str(len(x)))
#print("h length = " + str(len(h)))

# Timing variables for testing
reps = 5
resample1_times = [None] * reps
resample2_times = [None] * reps
slow_resample_times = [None] * reps

# Call the resample function reps times to get average times
for n in range(reps):
    y_fast, resample1_times[n] = resample1(x, h, U, D)
    y_fast2, resample2_times[n] = resample2(x, h, U, D)
    y_slow, slow_resample_times[n] = slow_resample(x, h, U, D)
    
print("Average time for resample1: " + str(list_avg(resample1_times) * 1000) + " ms")
print("Average time for resample2: " + str(list_avg(resample2_times) * 1000) + " ms")
print("Average time for slow_resample: " + str(list_avg(slow_resample_times) * 1000) + " ms")
print("")

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