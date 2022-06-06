import time


# 反向互补
def reverse_complement_all(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]


# 简化反向互补
def reverse_complement_limit(seq):
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


# kmers生成器
def make_kmers(seq, k):
    for i in range(len(seq) - k + 1): yield seq[i:i + k]


# 补齐生成器
def forward(seq):
    for x in 'ACGT': yield seq[1:] + x


def reverse(seq):
    for x in 'TGCA': yield x + seq[:-1]


def make_assemble_hashdict(fileslist, merSize, kmer_dict, ref_dict):
    kmer_count = 0
    for file in fileslist:
        infile = open(file, 'r', encoding='utf-8', errors='ignore')
        infile.readline()
        for line in infile:
            temp_str = []
            while line and line[0] != '>':
                temp_str.append(line)
                line = infile.readline()
            refseq_r = ''.join(filter(str.isalnum, ''.join(temp_str).upper()))
            for refseq in [refseq_r, reverse_complement_limit(refseq_r)]:
                kmer_count += len(refseq) - merSize + 1
                for j in range(0, len(refseq) - merSize + 1):
                    temp_list, temp_len, kmer = [], [0], refseq[j:j + merSize]
                    if kmer in kmer_dict:
                        temp_list = kmer_dict[kmer]
                        temp_list[0] += 1
                    else:
                        temp_list = [1, temp_len]
                        kmer_dict[kmer] = temp_list
                    if kmer in ref_dict:
                        temp_len[0] = 1 + ref_dict[kmer][1] / ref_dict[kmer][0] if ref_dict[kmer][1] < 0 else \
                        ref_dict[kmer][1] / ref_dict[kmer][0]
        infile.close()
    return kmer_count


def get_dis(_pos, new_pos, weight=0.5):
    return 1 - abs(_pos - new_pos) if (_pos and new_pos) else weight


def get_contig_forward_v5(_dict, seed, iteration=1000, weight=0.5):
    temp_list, reads_list, stack_list = [seed], [seed], []
    best_kmc, cur_kmc, best_pos, cur_pos, best_seq, cur_seq = [], [], [], [], [], []
    _pos, node_distance = 0, 0
    while True:
        node = sorted(
            [(i, _dict[i][1][0], _dict[i][0] ** get_dis(_pos, _dict[i][1][0], weight)) for i in forward(temp_list[-1])
             if i in _dict], key=lambda _: _[2], reverse=True)
        while node:
            if node[0][0] in temp_list:
                node.pop(0)
            else:
                break
        if not node:
            iteration -= 1
            if sum(cur_kmc) > sum(
                best_kmc): best_kmc, best_seq, best_pos = cur_kmc.copy(), cur_seq.copy(), cur_pos.copy()
            [(temp_list.pop(), cur_kmc.pop(), cur_seq.pop(), cur_pos.pop()) for _ in range(node_distance)]
            if stack_list:
                node, node_distance, _pos = stack_list.pop()
            else:
                break
        if len(node) >= 2:
            stack_list.append((node[1:], node_distance, _pos))
            node_distance = 0
        if node[0][1] > 0: _pos = node[0][1]
        temp_list.append(node[0][0])
        reads_list.append(node[0][0])
        cur_pos.append(node[0][1])
        cur_kmc.append(node[0][2])
        cur_seq.append(node[0][0][-1])
        node_distance += 1
        if not iteration: break
    return best_seq, reads_list, best_kmc, best_pos


def get_contig_v5(_dict, seed, iteration=1000, weight=0.5):
    right, reads_list_1, k_1, pos_1 = get_contig_forward_v5(_dict, seed, iteration, weight)
    left, reads_list_2, k_2, pos_2 = get_contig_forward_v5(_dict, reverse_complement_limit(seed), iteration, weight)
    _pos = [x for x in pos_1 + pos_2 if x > 0]
    min_pos, max_pos = 0, 1
    if _pos: min_pos, max_pos = min(_pos), max(_pos)
    return reverse_complement_limit(''.join(left)) + seed + ''.join(
        right), min_pos, max_pos, reads_list_1 + reads_list_2


def get_scaffold(_dict, seed_list, limit, ref_length, iteration=1000, weight=0.5):
    mid_seed, min_dis = seed_list[0], 1
    # for seed in seed_list:
    #     if abs(seed[2] - 0.5) < min_dis:
    #         mid_seed, min_dis = seed, abs(seed[2] - 0.5)

    contig, min_pos, max_pos, read_list = get_contig_v5(_dict, mid_seed[0], iteration, weight)

    for i in read_list:
        if i in _dict[i]: del _dict[i]
    contigs = [(contig, min_pos, max_pos)]
    new_seed = []

    for x in seed_list:
        if x[2] / x[1] < min_pos and x[1] > limit: new_seed = x
    while new_seed:
        contig, _min_pos, _max_pos, read_list = get_contig_v5(_dict, new_seed[0], iteration, weight)
        for i in read_list:
            if i in _dict[i]: del _dict[i]
        if new_seed[2] / new_seed[1] < _min_pos or new_seed[2] / new_seed[1] > _max_pos: break
        # print(_min_pos, _max_pos,contig)
        if _max_pos < contigs[0][2]:
            contigs.insert(0, (contig, _min_pos, _max_pos))
        else:
            break
        new_seed = []
        for x in seed_list:
            if x[2] / x[1] < _min_pos and x[1] > limit: new_seed = x
    new_seed = []
    for x in seed_list:
        if x[2] / x[1] > max_pos and x[1] > limit:
            new_seed = x
    while new_seed:
        contig, _min_pos, _max_pos, read_list = get_contig_v5(_dict, new_seed[0], iteration, weight)
        for i in read_list:
            if i in _dict[i]: del _dict[i]
        if new_seed[2] / new_seed[1] < _min_pos or new_seed[2] / new_seed[1] > _max_pos: break
        # print(_min_pos, _max_pos, contig)
        if _min_pos > contigs[-1][1]:
            contigs.append((contig, _min_pos, _max_pos))
        else:
            break
        new_seed = []
        for x in seed_list:
            if x[2] / x[1] > _max_pos and x[1] > limit: new_seed = x
    scaffold = contigs[0][0]
    for x in range(1, len(contigs)):
        scaffold += max(2, int((contigs[x][1] - contigs[x - 1][2]) * ref_length)) * 'N' + contigs[x][0]
    return scaffold


def dynamic_limit_v1(_dict, smoother_level=4, smoother_size=64):
    count_list = [0] * smoother_size
    for x in _dict:
        if _dict[x][0] <= smoother_size: count_list[_dict[x][0] - 1] += 1
    # 平滑器
    for i in range(smoother_level):
        for x in range(1, smoother_size - smoother_level + i):
            if count_list[x + smoother_level - 1 - i] < count_list[x - 1] and count_list[x] < count_list[
                x + smoother_level - i]:
                return x + 1
    return 2


def dynamic_limit_v2(_dict, ref_length, N, ref_rate=1.5, list_size=256):
    count_list = [0] * list_size
    for x in _dict:
        if _dict[x][0] <= list_size: count_list[_dict[x][0] - 1] += 1
    # 计算器
    F0, sum_f, sum_k = len(_dict), 0, 0
    for x in range(list_size):
        sum_f += count_list[x]
        sum_k += count_list[x] * (x + 1)
        if (F0 - sum_f) / 2 < ref_length * ref_rate: return max(2, x), round((N - sum_k) / ref_length / ref_rate / 2, 2)
    return 2, round((N - sum_k) / ref_length / ref_rate / 2, 2)


def dynamic_weight(_dict):
    temp_count = 0
    for i in _dict:
        if _dict[i][1][0] > 0: temp_count += 1
    return 1 - temp_count / len(_dict)


if __name__ == "__main__":
    limit = -1
    k_value = 20
    value = 1514
    ref_dict, filtered_dict = {}, {}

    from miner import gethashdict_v4

    gethashdict_v4(["./refs/Apiale/7128.fasta"], k_value, ref_dict, True, False, True)
    kmer_sum = make_assemble_hashdict(['./results_AP_35/7128.fasta'], k_value, filtered_dict, ref_dict)

    for i in ref_dict:
        if i not in filtered_dict:
            ref_dict[i][0] = 0

    if limit < 0:
        limit, deepth = dynamic_limit_v2(filtered_dict, value, kmer_sum, 2)
        if deepth > 20: limit = max(limit, dynamic_limit_v1(filtered_dict, 8, 64))
    print('\n', limit, deepth)

    if limit > 0:
        _filter = [x for x in filtered_dict if filtered_dict[x][0] <= limit]
        for x in _filter:
            if filtered_dict[x][1][0] == 0:
                del filtered_dict[x]
        _filter = []
    cur_weight = dynamic_weight(filtered_dict) if filtered_dict else 0.5
    print(cur_weight)

    seed_list = [(x, ref_dict[x][0], ref_dict[x][1]) for x in ref_dict if ref_dict[x][0] > 0 and ref_dict[x][1] > 0]
    list.sort(seed_list, key=lambda x: x[1], reverse=True)

    t1 = time.time()

    last_seed = seed_list[0][0]
    contig = get_contig_v5(filtered_dict, last_seed, weight=cur_weight)[0]
    # if len(contig) < value:
    #     for seed in seed_list[1:]:
    #         if last_seed[0][1:] == seed[0][:-1] or last_seed[0][2:] == seed[0][:-2]:
    #             last_seed = seed
    #             continue
    #         last_seed = seed
    #         contig = get_contig_v5(filtered_dict, seed[0])[0]
    #         if len(contig) > value:
    #             break

    # contig = get_scaffold(filtered_dict, seed_list, limit, value)
    print("\n>scaffold_" + str(len(contig)))
    print(contig)

    t2 = time.time()

    print("time:", t2 - t1)
