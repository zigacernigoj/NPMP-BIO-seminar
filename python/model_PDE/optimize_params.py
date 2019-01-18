import numpy as np
import time
import math

def get_sync_index(A_full, treshold, num_of_elements):

    start_time_getsync = time.time()

    # A_full[vrstice, stolpci] = A_full[casovni intervali, celice]
    A_diff = [np.max(A_full[i, :]) - np.min(A_full[i, :]) for i in range(A_full.shape[0])]
    A_lower = list(map(lambda x: x <= treshold, A_diff))

    trueArray = [True for _ in range(num_of_elements)]
    synced = [(i, i + num_of_elements) for i in range(len(A_lower)) if A_lower[i:i + num_of_elements] == trueArray]

    # print("razlike", A_diff)
    # print("is lower", A_lower)
    # print("synced", synced[0][0])

    # print("--- %s seconds for get sync ---" % (time.time() - start_time_getsync))

    if len(synced) > 0:
        sync_index = synced[0][0]
        min_for_cell = np.min(A_full[sync_index:, 0])
        max_for_cell = np.max(A_full[sync_index:, 0])

        if (max_for_cell - min_for_cell) < 150:
            return math.inf
        else:
            return synced[0][0]
    else:
        return math.inf


# def getT(A_lower, num_of_elements):
#     i = 0
#     counter = 0
#     while i < len(A_lower):
#         if A_lower[i]:
#             counter += 1
#         else:
#             counter = 0
#         if counter == num_of_elements:
#             print("--- %s seconds for get sync ---" % (time.time() - start_time_getsync))
#             return i - num_of_elements + 1
#         i += 1
#
#     print("--- %s seconds for get sync ---" % (time.time() - start_time_getsync))
#     return -1






