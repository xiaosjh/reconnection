#include <stdio.h>

int main() {
    double va = 0.5, vb = 0.3;
    
    // 计算第一次相遇的时间
    double first_meeting_time = 2 * 2 / (va + vb);
    printf("第一次相遇的时间为 %.2f h\n", first_meeting_time);
    
    // 计算第二次相遇的时间
    double second_meeting_time = 4 * 2 / (va + vb);
    printf("第二次相遇的时间为 %.2f h\n", second_meeting_time);
    
    return 0;
}