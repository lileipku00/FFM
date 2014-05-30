import matplotlib.pyplot as plt

# 70 station-receiver pairs
x = [1, 2, 3, 4]
y = [62.4, 36.1, 34.7, 33.3]
plt.plot(x, y, 'b', lw=3, label='70')

# 172 station-receiver pairs
x = [1, 2, 3, 4]
y = [128, 71.4, 68.6, 62.5]
plt.plot(x, y, 'r', lw=3, label='172')

# 869 station-receiver pairs
x = [1, 2, 3, 4]
y = [0, 5.*60 + 29.6, 4.*60 + 55.5, 4.*60+35.4]
plt.plot(x, y, 'k', lw=3, label='869')


plt.xlabel('#Processes', size='x-large', weight='bold')
plt.ylabel('Time (sec)', size='x-large', weight='bold')

plt.xticks([1,2,3,4], size='x-large', weight='bold')
plt.yticks(size='x-large', weight='bold')

plt.legend()
plt.show()
