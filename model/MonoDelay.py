import numpy as np

def delayBlock(input_block, state_block):
    output_block = np.concatenate((state_block, input_block[:-len(state_block)]))
    state_block = input_block[-len(state_block):]
    return output_block, state_block

def generateBlock(block_size):
    out_block = np.random.rand(block_size)
    for n in range(len(out_block)):
        out_block[n] = round(out_block[n] * 100)
    return out_block

block = generateBlock(10)
state = generateBlock(4)
out_block, out_state = delayBlock(block, state)


for n in range(len(block)):
    print("input block[" + str(n) + "] = " + str(block[n]))

print("")

for n in range(len(state)):
    print("input state[" + str(n) + "] = " + str(state[n]))

print("")

for n in range(len(block)):
    print("output block[" + str(n) + "] = " + str(out_block[n]))

print("")

for n in range(len(state)):
    print("output state[" + str(n) + "] = " + str(out_state[n]))

print("")
