#include <stdio.h>
#include <string.h>
#include <stdlib.h>


/***
 utils.c 用来实现Easy353的部分算法或者函数 用于加快速度
        最终以动态链接库的方式让python调用
***/
int max(int a, int b) {
    return a > b ? a : b;
}

// 接受两个字符串s1和s2作为输入
// 计算两个字符串之间的最长公共子串的长度 longest common substring
// 使用动态规划算法求解，并通过滚动数组优化空间复杂度
// min_overlap_len 和 max_overlap_len是对最长公共子串的长度进行限制
int optimized_lcs_len(char *s1, char *s2,int min_overlap_len,int max_overlap_len) {
    int max_len = 0;
    int len_s1 = strlen(s1);
    int len_s2 = strlen(s2);
    int *prev_col = (int *)calloc(len_s2 + 1, sizeof(int));
    int *curr_col = (int *)calloc(len_s2 + 1, sizeof(int));
    int i, j;

    for (i = 1; i <= len_s1; i++) {
        for (j = 1; j <= len_s2; j++) {
            if (s1[i - 1] == s2[j - 1]) {
                curr_col[j] = prev_col[j - 1] + 1;
                max_len = max(max_len, curr_col[j]);
                if(max_len >= max_overlap_len){
                    free(prev_col);
                    free(curr_col);
                    return max_len;
                }
            } else {
                curr_col[j] = 0;
            }
        }
        int *temp = prev_col;
        prev_col = curr_col;
        curr_col = temp;
    }

    free(prev_col);
    free(curr_col);
    return max_len >= min_overlap_len ? max_len : 0;
}

int main(){
    char s1[] = "AATGTTGTTCAAGATGAGGACAAGCTTGTCACTTCAAACACTGATTGGATGCATAAATACAAAGGCTCCAGTAAGCTCCTGTTACAACCTAGGACCGCTGATCAGGTTTCTCAGATTCTTAAATATTGTAACTCCAGAAACTTGGCTGTTGTTCCCCAAGGTGGAAATACTGGTCTTGTAGGAGGAAGTGTGCCTGTCTTTGATGAAGTGATTGTTAGTCTTAGTTCTATGAATAAGATCATATCTTTCGACAAGGTTAGTGGAATATTGGTATGTGAAGCTGGGTGCATATTGGAAAATATAATGTCATTCCTGGACAATGAAGGATTTATTATGCCACTAGACTTAGGTGCAAAAGGGAGTTGTCAGATTGGTGGAAATGTTTCAACTAATGCTGGTGGTTTGCGTCTTGTTCGCTATGGATCGCTTCATGGAAGTGTACTTGGTGTTGAAGCTGTTCTGGCAAATGGTACTGTACTTGACATGCTTAAGACATTGCGCAAAGATAATACTGGCTATGATTTGAAACATCTATTTATAGGAAGTGAAGGTTCCTTGGGAATTGTTACGAAGGTTTCAATACTTACCCCACCAAAGTTATCTTCAGTGAATGTGGCTTTTCTTGCTTGTAAAGACTATAGCAGCTGCCAGAAATTGCTACAGGAGGCAAAAGGGAAACTTGGGGAGATTTTATCTGCATTTGAATTTCTGGATGTCCAGTCTATGAATTTGGTTTTAAATCACATGGAAGGTGCACGAAATCCATTGCCGTCATTGCATAACTTTTATGTTTTGATTGAGACAACAGGCAGTGATGAATCTTCTGACAAGCAAAAGCTTGAAGCATTTCTACTTGGCTCCATGGAGAATGAATTGATATCCGATGGTGTTCTTGCACAAGACATAAACCAAGCATCGTCTTTTTGGCTTCTACGTGAGGGT";
    char s2[] = "ATAAATAGATGTTTCAAATCATAGCCAGTATTATCTTTGCGCAATGTCTTAAGCATGTCAAGTACAGTACCATTTGCCAGAACAGCTTCAACACCAAGTA";
    int result = optimized_lcs_len(s1, s2,1,8);
    printf("Length of Longest Common Substring: %d\n", result);
    return 0;
}