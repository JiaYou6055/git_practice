//LDPC Decoder EEP/UEP QC Mode  flag -1 with Random encode
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define CHK_MODE  1 //check node operation mode
                    // 0 = original delta= ln([1+exp^(-|L1+L2|)]/[1+exp^(-|L1-L2|)])
                    // 1 = table approximation
                    // 2 = min-sum delta=0
                    // 3 = modified min-sum

#define TABLE_BIT 2 //table approximation bit
                    // 0 = 2 bits  4 區段
                    // 1 = 3 bits  8 區段
                    // 2 = 4 bits  16區段

#define MODIFIED_A      0.9     //modified min-sum factor
#define ITERATION_MAX   100

#define STOP_CRITERION      1   //stop criterion 停止條件
                                // 0 = total error set
                                // 1 = UEP level-1 error set    (UEP case)
                                // 2 = UEP level-2 error set    (UEP case)

#define ERRONEOUS_SET_MAX   50  //累積 erroneous codeword 組數上限

#define LEVEL_1_BLK     7       //UEP Level-1
#define LEVEL_2_BLK     7       //UEP Level-2

#define ROW         3992
#define COLUMN      6986
#define I_MATRIX    499
#define INFO_BIT    2994            //encode part G的row
#define COL_MAX     1535            //encode part G的COL最多有#個1
#define CODE_RATE   0.4286
#define SNRdB       1.0          //SNR Eb/N0 (dB)


void Parameter_Display(char* input_file);
char* Get_Time(char* buffer,int mod);  // mod1 time display    mod2 txt output
void normal(double*n1 ,double*n2 ,double sigma);
double Ranq1();
double CHK_table(double x);
double CHK(double L1 , double L2);
void Encoder(int* codeword, int Info_length, int** G_matrix, int G_row, int G_col);

int** int_Matrix(int row, int col, int initial_value);          //Dynamic 宣告 int matrix
double** double_Matrix(int row, int col, double initial_value); //Dynamic 宣告 double matrix
int* int_Array(int ary_size,int initial_value);                 //Dynamic 宣告 int array
double* double_Array(int ary_size, double initial_value);       //Dynamic 宣告 double array

unsigned long long SEED = 104064521;    // SEED must be an unsigned integer smaller than 4101842887655102017.
unsigned long long RANV;
int RANI = 0;


//char qc_H_matrix[]="qc_uep_matrix_3786x7572_I_631_R_0.5007.txt";    // 6  6
//char qc_G_matrix[]="[encode]G_3790x7572_col_index_MAX1983.txt";

char qc_H_matrix[]="qc_uep_matrix_3992x6986_I_499_R_0.4286.txt";      // 7  7
char qc_G_matrix[]="[encode]G_2994x6986_col_index_MAX1535.txt";


int main()
{
    clock_t start = clock();    //timer
//------------------------------------------------------------------
//                variable initial
//------------------------------------------------------------------
    Parameter_Display(qc_H_matrix);

    FILE *fp=NULL;
    int i,j,k;
    int** qc_parity_matrix=int_Matrix(ROW/I_MATRIX,COLUMN/I_MATRIX,0);  //H( contracted version )
    int** parity_edge_row=int_Matrix(ROW, COLUMN/I_MATRIX+1 ,-1);       //parity check matrix ( +1是為了尾端flag:-1 )
    int** parity_edge_column=int_Matrix(ROW/I_MATRIX+1, COLUMN, -1);    //parity check matrix ( +1是為了尾端flag:-1)
    int mapping_index;
    int** col_row_covert=int_Matrix(ROW/I_MATRIX, COLUMN, -1);          //mapping
    int** row_col_covert=int_Matrix(ROW, COLUMN/I_MATRIX, -1);          //mapping

    int** g_matrix_column=int_Matrix(COLUMN, COL_MAX+1, -1);            //encode generator matrix (COL_MAX+1 多算-1 flag)
    int codeword[COLUMN]={0};
    //double codeword[COLUMN]={0};

    double* rx_codeword=double_Array(COLUMN,0);     //mapping & received C y=x+n
    double n1,n2;                                   //awgn noise
    double snr=pow(10,(double)SNRdB/10);            //SNR Eb/N0
    double sigma=sqrt(1/(2*CODE_RATE*snr));         //standard deviation
    double* initial_Ll=double_Array(COLUMN,0);      //initialization Ll=(4Es/N0)*yl
    //message passing part
    int count_iteration=0;
    //step 1
    double L1_tmp=0,L2_tmp=0;
    double** MP_chk=double_Matrix(ROW,COLUMN/I_MATRIX,0);
    int flag_chk=0;
    //step 2
    double L_tmp=0;
    double** MP_var=double_Matrix(ROW/I_MATRIX,COLUMN,0);
    double* total_llr=double_Array(COLUMN,0);
    //Hard decision part
    int* decode_codeword=int_Array(COLUMN,0);
    //syndrome part
    int parity_eqtmp=0;
    int flag_syndrome=0;
    //BER part
    unsigned int count_codeword_set=0;          //count total sended codeword  組數
    unsigned int count_total_error_set=0;       //count num total error codeword
    unsigned int count_stop_error_set=0;        //count num stop error codeword  (停止條件)
    unsigned int count_total_error_bits=0;      //count error bit 數 (累積)
    unsigned int count_level_error_bits=0;      //count level error bit 數 (累積)
    int flag_allcorrect=0;                      //decode all correct flag
    int flag_total_allcorrect=0;
    double ber_uep=0;                               //bit error rate
    double ber_toatal=0;                        //total bit error rate
//------------------------------------------------------------------
//          file open & input H & G
//------------------------------------------------------------------
    if( (fp=fopen(qc_H_matrix,"r")) ==NULL){
        printf("cannot open QC Mode file");
        exit(1);
    }
    while(!feof(fp)){
        for(i=0;i<ROW/I_MATRIX;i++){
            for(j=0;j<COLUMN/I_MATRIX;j++){
                fscanf(fp,"%d",&qc_parity_matrix[i][j]);
            }
        }
        break;
    }
    fclose(fp);

    if( (fp=fopen(qc_G_matrix,"r")) ==NULL){
        printf("cannot open Generator matrix file");
        exit(1);
    }
    while(!feof(fp)){
        for(i=0;i<COLUMN;i++){
            for(j=0;j<COL_MAX+1;j++){
                fscanf(fp,"%d",&g_matrix_column[i][j]);
                if(g_matrix_column[i][j]==-1)
                    break;
            }
        }
        break;
    }
    fclose(fp);

//-------------------------------------------------
//          generate QC form (若有I_0要注意)
//-------------------------------------------------
    int shift_row_tmp=0;
    int shift_col_tmp=0;
    int* row_edge_pointer   = int_Array(ROW/I_MATRIX,0);        //避免讀到0後沒有連續寫入row_edge裡
    int* column_edge_pointer= int_Array(COLUMN/I_MATRIX,0);     //紀錄此row/column現在放到的位置

    for(i=0;i<ROW/I_MATRIX;i++){            //  qc row
        for(j=0;j<COLUMN/I_MATRIX;j++){     //  qc column

            if(qc_parity_matrix[i][j]==0)
                continue;
            shift_row_tmp=qc_parity_matrix[i][j];     //Ie  e: 0 <= i < I_MATRIX 但0要特別處理,因為0目前被當成O(all-zero)

            for(k=0;k<I_MATRIX ;k++){
                //row_edge [row][?]
                //               整體偏移量          存入輔助指示           cyclic 位移           整體偏移量
                parity_edge_row[ I_MATRIX*i+ k ][ row_edge_pointer[i] ]= shift_row_tmp%I_MATRIX + I_MATRIX*j;

                //column_edge [?][column]
                parity_edge_column[ column_edge_pointer[j] ][ I_MATRIX*j+ shift_row_tmp%I_MATRIX ]= shift_col_tmp%I_MATRIX + I_MATRIX*i;

                shift_row_tmp++;        //qc 特性:右移1
                shift_col_tmp++;        //qc 特性:右移1
            }
            shift_row_tmp=0;            //reset
            shift_col_tmp=0;            //reset
            row_edge_pointer[i]++;      //紀錄每塊大row現在存到哪裡
            column_edge_pointer[j]++;   //紀錄每塊大column現在存到哪裡
        }
    }

    free(qc_parity_matrix[0]);
    free(qc_parity_matrix);
    free(row_edge_pointer);
    free(column_edge_pointer);

//-------------------------------------------------
//      mapping table col row
//                          col[][A]=B row[B][]=A
//-------------------------------------------------
    //column to row covert mapping                              找col方向看下來1的位置對應到row角度在哪 row[?][?]
    for(j=0;j<COLUMN;j++)                                       //control column index
    {
        for(i=0;i<(ROW/I_MATRIX+1);i++)                         //control Gamma index
        {
            mapping_index=parity_edge_column[i][j];             //col[i][A]=B  : 表第i個1出現在第B個row
            if(mapping_index==-1)                               //parity_edge_column 讀到flag -1
                break;
            for(k=0;k<(COLUMN/I_MATRIX+1);k++)
            {                                                   //find row[B][?]=A :找出第B個row出現第?個1時在第A個col
                if(parity_edge_row[mapping_index][k]==j)
                    col_row_covert[i][j]=k;                     //記錄下col[i][A]轉換到row是第?個1 (row[B][?]=A)
            }                                                   //row[ col[i][A] ][ col_row_covert[i][A] ]
        }
    }
    //row to column covert mapping                              找row方向往右看1的位置對應到col角度在哪 col[?][?]
    for(i=0;i<ROW;i++)                                          //control ROW index
    {
        for(j=0;j<(COLUMN/I_MATRIX+1);j++)                      //control RHO index
        {
            mapping_index=parity_edge_row[i][j];                //row[B][?]=A  : 表第i個1出現在第A個column
            if(mapping_index==-1)                               //parity_edge_row 讀到flag -1
                break;
            for(k=0;k<(ROW/I_MATRIX+1);k++)
            {                                                   //find col[?][A]=B :找出第A個column出現第?個1時在第B個ROW
                if(parity_edge_column[k][mapping_index]==i)
                    row_col_covert[i][j]=k;
            }
        }
    }

//------------------------------------------------------------------
//          generate codeword
//------------------------------------------------------------------
    srand(SEED);    //為了在相同的seed下值會一樣

RESEND_CODEWORD:
    Encoder(codeword,INFO_BIT,g_matrix_column,COLUMN,COL_MAX+1);    //random info encode codeword

//------------------------------------------------------------------
//          mapping & addition noise   y=x+n    0 --> 1   1 ---> -1
//------------------------------------------------------------------
//RESEND_CODEWORD:
    for(i=0;i<COLUMN;i=i+2){
        normal(&n1,&n2,sigma);
        rx_codeword[i]    = (codeword[i]  == 0 ? 1: -1)+n1;
        rx_codeword[i+1]  = (codeword[i+1]== 0 ? 1: -1)+n2;
//        printf("%lf\n%lf\n",rx_codeword[i],rx_codeword[i+1]);
    }

//------------------------------------------------------------------
//          Initialization LLR
//------------------------------------------------------------------
    for(i=0;i<COLUMN;i++){
        initial_Ll[i]=4*CODE_RATE*snr*rx_codeword[i];
    }
//------------------------------------------------------------------
//          Message passing
//------------------------------------------------------------------

    count_iteration=0;                  //message passing iteration count reset
    while(count_iteration<ITERATION_MAX)
    {
        //-------------------------------------------------
        //      step1 check node
        //-------------------------------------------------
        for(k=0;k<ROW;k++){                         //control row index
            for(i=0;i<COLUMN/I_MATRIX;i++){         //control rho index

                if(parity_edge_row[k][i]==-1)       //[QC] 碰到edge flag -1
                    break;

                for(j=0;j<COLUMN/I_MATRIX;j++){

                    if(parity_edge_row[k][j]==-1)                       //[QC] 碰到edge table底 flag -1 不用做了
                        break;

                    if((parity_edge_row[k][i])==(parity_edge_row[k][j]))    //避開自己
                        continue;

                    if(flag_chk==0){                                     //把第一筆當作L1,flag_chk=-1 標記為已有第一個輸入

                        if (count_iteration==0){                        //count_iteration=0 第一次做遞迴
                            L1_tmp=initial_Ll[parity_edge_row[k][j]];   //第一次遞迴 用LLR

                        }else{
                            L1_tmp=MP_var[ row_col_covert[k][j] ][ parity_edge_row[k][j] ];

                        }
                        flag_chk=-1;
                        continue;
                    }

                    if (count_iteration==0){                            //count_iteration=0 第一次做遞迴
                        L2_tmp=initial_Ll[parity_edge_row[k][j]];       //第一次遞迴 用LLR

                    }else{
                        L2_tmp=MP_var[ row_col_covert[k][j] ][ parity_edge_row[k][j] ];

                    }

                    L1_tmp=CHK(L1_tmp,L2_tmp);
                }
                flag_chk=0;                 //reset CHK flag
                MP_chk[k][i]=L1_tmp;
            }
        }
        //-------------------------------------------------
        //      step2 variable node & step 3 total llr
        //               col[][A]=B row[B][]=A
        //-------------------------------------------------
        for(i=0;i<COLUMN;i++){
            for(j=0;j<ROW/I_MATRIX;j++){

                if(parity_edge_column[j][i]==-1)        //[QC] 碰到edge flag -1 不用做了
                    break;

                for(k=0;k<ROW/I_MATRIX;k++){
                    if(parity_edge_column[k][i]==-1)    //[QC] 碰到edge table底 flag -1 不用做了
                        break;

                    if(parity_edge_column[k][i]==parity_edge_column[j][i])  // avoid 指到自己
                    {
                        L_tmp += initial_Ll[i];                             // 加上initial Ll 因為一定會碰到自己 把初值放在這裡加
                        total_llr[i] += (initial_Ll[i] + MP_chk[parity_edge_column[k][i]][col_row_covert[k][i]]);
                        continue;
                    }
                    L_tmp += MP_chk[parity_edge_column[k][i]][col_row_covert[k][i]];
                    total_llr[i] += MP_chk[parity_edge_column[k][i]][col_row_covert[k][i]];
                }
                MP_var[j][i]=L_tmp;
                L_tmp=0;            //reset L_tmp  讓下一次variable node 運算能正常
            }
        }
        //-------------------------------------------------
        //  Hard decision  1 if llr >= 0 ; 0 if llr < 0
        //-------------------------------------------------
        for(i=0;i<COLUMN;i++){
            if(total_llr[i]>=0){
                decode_codeword[i]=0;
            }else{
                decode_codeword[i]=1;
            }
            total_llr[i]=0;         //llr reset 避免在step 2&3 保留前一次iteration的值
        }
        //-------------------------------------------------
        //  Syndrome  Hx=0
        //-------------------------------------------------
        flag_syndrome=0;                                //syndrome 錯誤flag  0-->syndrome=0  -1 --> syndrome有1
        for(j=0;j<ROW;j++)
        {
            for(i=0;i<COLUMN/I_MATRIX;i++)
            {
                if(i==0)                                //[[ warning ]]每個row一定會有一個1 就不檢查是否碰到flag -1
                {
                    parity_eqtmp=decode_codeword[parity_edge_row[j][i]];
                    continue;
                }

                if(parity_edge_row[j][i]==-1)           //碰到flag -1 :此row已經到edge底了 沒有1
                    break;

                parity_eqtmp^=decode_codeword[parity_edge_row[j][i]];        //XOR : sum of parity check equation=0
            }
            if(parity_eqtmp==1)         //syndrome 有1 --> error
            {
                flag_syndrome=-1;       //syndrome 有錯 flag
                break;
            }
        }
        if (flag_syndrome==0)           //decode 正確-->結束迴圈
            break;
        else
            count_iteration++;
    }
    count_codeword_set++;               //count total sended codeword : 數完成一組codeword傳送

    //-------------------------------------------------
    //          Error bit count
    //-------------------------------------------------
    flag_allcorrect=0;          //flag reset            error --> -1   all correct --> 0
    flag_total_allcorrect=0;    //flag total error set  error --> -1   all correct --> 0

    for(i=0;i<COLUMN;i++)       //count error bits
    {
        if(codeword[i]!=decode_codeword[i]){    //total error set general case

            flag_total_allcorrect=-1;           //只要codeword有錯就立flag
            count_total_error_bits++;

            if(STOP_CRITERION==0)
                flag_allcorrect=-1;             //codeword有錯 且 滿足區段條件

            if( (STOP_CRITERION==1) && (i<LEVEL_1_BLK*I_MATRIX) ){                              // UEP level-1 error set
                count_level_error_bits++;
                flag_allcorrect=-1;             //codeword有錯 且 錯在某區段
            }
            if( (STOP_CRITERION==2) && ( LEVEL_1_BLK*I_MATRIX<=i && i<(LEVEL_1_BLK+LEVEL_2_BLK)*I_MATRIX) ){  // UEP level-2 error set
                count_level_error_bits++;
                flag_allcorrect=-1;             //codeword有錯 且 錯在某區段
            }
        }
    }

    if(flag_total_allcorrect==-1)
        count_total_error_set++;        //total error set

    if(flag_allcorrect==-1){            //STOP CRITERION count 累積error的codeword數

        count_stop_error_set++;

        if(STOP_CRITERION==0){          //total error set general case
            printf("[COUNT] Total error set = %u\n",count_stop_error_set);

        }else if(STOP_CRITERION==1 || STOP_CRITERION==2){                     // UEP level-1 and level-2 error set
            printf("[COUNT] Level-%d error set = %u\n",(STOP_CRITERION==1?1:2),count_stop_error_set);
            if( (count_stop_error_set%10) ==0)
                printf("[Info] Level-%d BER = %0.5e\n",(STOP_CRITERION==1?1:2),
                    (double)count_level_error_bits /( (double)( STOP_CRITERION==1?LEVEL_1_BLK:LEVEL_2_BLK )* I_MATRIX* count_codeword_set) );

        }else{
            printf("[ERROR] STOP_CRITERION value is fault!\n");
            exit(1);
        }

//        printf("[UEP Count] Level-%d error bits = %u\n",(STOP_CRITERION==1?1:2),count_level_error_bits);
//        printf("[Count]  Total  error bits     = %u\n",count_total_error_bits);

        if( (count_stop_error_set%10) ==0)
            printf("[Info] Total   BER = %0.5e\n",(double)count_total_error_bits/((double)count_codeword_set*(double)COLUMN));
    }

    if(count_stop_error_set<ERRONEOUS_SET_MAX)     //count_stop_error_set codeword error 組數計數
    {
        goto RESEND_CODEWORD;
    }
    //-------------------------------------------------
    //          bit error rate
    //-------------------------------------------------
    ber_toatal=(double)count_total_error_bits/((double)count_codeword_set*(double)COLUMN);
    ber_uep=(double)count_level_error_bits/( (double)( ( STOP_CRITERION==1 ? LEVEL_1_BLK : LEVEL_2_BLK ) * I_MATRIX * count_codeword_set ) );

    //-------------------------------------------------
    printf("\n----------[Result]----------\n\n");
    if(CHK_MODE==0){
        printf("[MODE=%d] original delta\n",CHK_MODE);
    }
    else if(CHK_MODE==1){
        printf("[MODE=%d] table approximation\n",CHK_MODE);
    }
    else if(CHK_MODE==2){
        printf("[MODE=%d] delta=0\n",CHK_MODE);
    }else if(CHK_MODE==3){
        printf("[MODE=%d] modified min-sum [a = %lf ]\n",CHK_MODE,MODIFIED_A);
    }

    printf("[SEED]          %llu\n",SEED);
    printf("[Iteration]     %d\n",ITERATION_MAX);
    printf("[Erroneous set] %d\n",ERRONEOUS_SET_MAX);
    printf("[Count]  Total  Tx_C  set      = %u\n",count_codeword_set);
    printf("[Count]  Total  error bits     = %u\n",count_total_error_bits);

    if(STOP_CRITERION==0){
    printf("SNRdB  = %lf\n",SNRdB);
    printf("BER    = %0.5e\n",ber_toatal);
    }
    else if(STOP_CRITERION==1 || STOP_CRITERION==2){
    printf("[Count] UEP Level-%d error bits = %u\n",(STOP_CRITERION==1?1:2),count_level_error_bits);
    printf("SNRdB  = %lf\n",SNRdB);
    printf("Total BER          = %0.5e\n",ber_toatal);
    printf("[UEP] Level-%d BER  = %0.5e\n",(STOP_CRITERION==1?1:2),ber_uep);
    }

    //timer
    clock_t end = clock();                              // 宣稱資料型態為 clock_t 的 end 變數
    int costTime = (float)(end - start)/CLK_TCK;        // 求算CPU執行時間
    printf("[Info.] Operation Time  %d(h) %d(m) %d(s)\n\n",costTime/3600,(costTime/60)%60, costTime%60);

    //-----------------------------------------------------------
    //          Save data
    //-----------------------------------------------------------
    char txtname[40];
    if( (fp=fopen(Get_Time(txtname,2),"w")) ==NULL){  //"\backup\ " Get_Time(txtname,2)
        printf("cannot generate save data file");
        exit(1);
    }
    fprintf(fp,"%d x %d\tI=%d\tR=%.5f\n",ROW,COLUMN,I_MATRIX,CODE_RATE);
    if(CHK_MODE==0){
        fprintf(fp,"[MODE=%d] original delta\n",CHK_MODE);
    }
    else if(CHK_MODE==1){
        fprintf(fp,"[MODE=%d] table approximation\n",CHK_MODE);
    }
    else if(CHK_MODE==2){
        fprintf(fp,"[MODE=%d] delta=0\n",CHK_MODE);
    }else if(CHK_MODE==3){
        fprintf(fp,"[MODE=%d] modified min-sum [a = %lf ]\n",CHK_MODE,MODIFIED_A);
    }
    fprintf(fp,"[SEED]\t%llu\n",SEED);
    fprintf(fp,"[Iteration]\t%d\n",ITERATION_MAX);
    fprintf(fp,"[Erroneous set]\t%d\n",ERRONEOUS_SET_MAX);
    fprintf(fp,"[Count]  Total  Tx_C  set\t= %u\n",count_codeword_set);
    fprintf(fp,"[Count]  Total  error bits\t= %u\n",count_total_error_bits);
    if(STOP_CRITERION==0){
    fprintf(fp,"SNRdB\t= %lf\n",SNRdB);
    fprintf(fp,"BER\t\t= %0.5e\n",ber_toatal);
    }
    else if(STOP_CRITERION==1 || STOP_CRITERION==2){
    fprintf(fp,"[Count] UEP Level-%d error bits\t= %u\n",(STOP_CRITERION==1?1:2),count_level_error_bits);
    fprintf(fp,"SNRdB  = %lf\n",SNRdB);
    fprintf(fp,"Total BER\t\t= %0.5e\n",ber_toatal);
    fprintf(fp,"[UEP] Level-%d BER\t= %0.5e\n",(STOP_CRITERION==1?1:2),ber_uep);
    }
    fprintf(fp,"[Info.] Operation Time  %d(h) %d(m) %d(s)\n\n",costTime/3600,(costTime/60)%60, costTime%60);
    fclose(fp);

    printf("\a");
    system("pause");

    //------------------------------------------------------------
    //          free
    //------------------------------------------------------------
    free(parity_edge_row[0]);
    free(parity_edge_row);
    free(parity_edge_column[0]);
    free(parity_edge_column);
    free(row_col_covert[0]);
    free(row_col_covert);
    free(col_row_covert[0]);
    free(col_row_covert);
    free(rx_codeword);
    free(initial_Ll);
    free(MP_chk[0]);
    free(MP_chk);
    return 0;
}

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
//
//                  Function
//
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
double CHK(double L1 ,double L2){       // initial_Ll[index 1] , initial[index 2]
    int a;
    double modified_factor=1;
    double r,delta;

    if(CHK_MODE==0){
        delta=log((1+exp(-1*fabs(L1+L2)))/(1+exp(-1*fabs(L1-L2)))); //ln([1+exp^(-|L1+L2|)]/[1+exp^(-|L1-L2|)])
    }
    else if(CHK_MODE==1){
        delta=CHK_table(L1+L2)-CHK_table(L1-L2);
    }
    else{
        delta=0;     //min-sum algorithm
    }

    if(CHK_MODE==3){                    //modified min-sum
        modified_factor=MODIFIED_A;
    }else{
        modified_factor=1;
    }

    a=(L1*L2)>0?1:-1;       //sgn(L1)*sgn(L2)

    if (fabs(L1)<fabs(L2)){  //min(|L1|,|L2|)
        r=modified_factor*a*fabs(L1)+delta;
    }
    else{
        r=modified_factor*a*fabs(L2)+delta;
    }
    return r;
}

void normal(double *n1 ,double *n2 ,double sigma){
    double x1,x2,s;
    do{
        x1=Ranq1();
        x2=Ranq1();
        x1=2*x1-1;
        x2=2*x2-1;
        s=x1*x1+x2*x2;
    }while(s>=1.0);
    *n1=sigma*x1*sqrt(-2*log(s)/s);
    *n2=sigma*x2*sqrt(-2*log(s)/s);
}

double Ranq1(){
    if ( RANI == 0 ){
        RANV = SEED ^ 4101842887655102017LL;
        RANV ^= RANV >> 21;
        RANV ^= RANV << 35;
        RANV ^= RANV >> 4;
        RANV = RANV * 2685821657736338717LL;
        RANI++;
    }
    RANV ^= RANV >> 21;
    RANV ^= RANV << 35;
    RANV ^= RANV >> 4;
    return RANV * 2685821657736338717LL * 5.42101086242752217E-20;
}

double CHK_table(double x){
    switch(TABLE_BIT){
        case 0:
            //4區間
            if(-0.56<x && x<0.56)
                return 0.57245;
            else if( (0.56<=x && x<1.5) || (-0.56>=x && x>-1.5) )
                return 0.3266;
            else if( (1.5 <=x && x<4.6) || (-1.5 >=x && x>-4.6) )
                return 0.1057;
            else if (4.6<=x || -4.6>=x)
                return 0;
            break;
        case 1:
            //8區間
            if(-0.19<x && x<0.19)
                return 0.6479;
            else if( (0.19<=x && x<0.43) || (-0.19>=x && x>-0.43) )
                return 0.5519;
            else if( (0.43<=x && x<0.7) || (-0.43>=x && x>-0.7))
                return 0.45215;
            else if( (0.7<=x && x<1.05) || (-0.7>=x && x>-1.05))
                return 0.35165;
            else if( (1.05<=x && x<1.5) || (-1.05>=x && x>-1.5))
                return 0.25075;
            else if( (1.5<=x && x<2.25) || (-1.5>=x && x>-2.25))
                return 0.1508;
            else if( (2.25<=x && x<4.6) || (-2.25>=x && x>-4.6))
                return 0.0551;
            else if(4.6<=x || -4.6>=x)
                return 0;
            break;
        case 2:
            //16區間
            if(-0.19<x && x<0.19)
                return 0.6479;
            else if( (0.19<=x && x<0.43) || (-0.19>=x && x>-0.43) )
                return 0.5519;
            else if( (0.43<=x && x<0.7) || (-0.43>=x && x>-0.7))
                return 0.45215;
            else if( (0.7<=x && x<1.05) || (-0.7>=x && x>-1.05))
                return 0.35165;
            else if( (1.05<=x && x<1.5) || (-1.05>=x && x>-1.5))
                return 0.25075;
            else if( (1.5<=x && x<2.25) || (-1.5>=x && x>-2.25))
                return 0.1508;
            else if( (2.25<=x && x<2.97) || (-2.25>=x && x>-2.97))
                return 0.07512;
            else if( (2.97<=x && x<3.97) || (-2.97>=x && x>-3.97))
                return 0.0376;
            else if( (3.97<=x && x<4.6) || (-3.97>=x && x>-4.6))
                return 0.01758;
            else if(4.6<=x || -4.6>=x)
                return 0;
            break;
    }
}
//Dynamic 宣告 int matrix
int** int_Matrix(int row, int col, int initial_value){
    int i=0,j=0;
    int **ary_r = (int**)calloc( row,sizeof(int*));
    if ( ary_r==NULL ){                         //避免要不到記憶體
        printf("[warning.]Allocate fault !\n");
        return NULL;
    }
    int *pData = (int*)calloc(row*col,sizeof(int));
    if ( pData==NULL ){                         //避免要不到記憶體
        printf("[warning.]Allocate fault !\n");
        free(ary_r);                            //失敗的話前面要的記憶體也要釋放
        return NULL;
    }
    for(i = 0; i < row ; i++, pData += col)
        ary_r[i] = pData;

    if(initial_value==0)
        return ary_r;
    else{
        for(i=0; i< row; i++)
            for(j=0 ; j< col ;j++)
                ary_r[i][j]=initial_value;
        return ary_r;
    }
}
//Dynamic 宣告 double matrix
double** double_Matrix(int row, int col, double initial_value){
    int i=0,j=0;

    double **ary_r = (double**)calloc( row,sizeof(double*));
    if ( ary_r==NULL ){                         //避免要不到記憶體
        printf("[warning.]Allocate fault !\n");
        return NULL;
    }
    double *pData = (double*)calloc(row*col,sizeof(double));
    if ( pData==NULL ){                         //避免要不到記憶體
        printf("[warning.]Allocate fault !\n");
        free(ary_r);                            //失敗的話前面要的記憶體也要釋放
        return NULL;
    }
    for(i = 0; i < row ; i++, pData += col)
        ary_r[i] = pData;

    if(initial_value==0)
        return ary_r;
    else{
        for(i=0; i< row; i++)
            for(j=0 ; j< col ;j++)
                ary_r[i][j]=initial_value;
        return ary_r;
    }
}
//Dynamic 宣告 int array
int* int_Array(int ary_size, int initial_value){
    int i=0;
    int* ary=calloc(ary_size,sizeof(int));
    if(ary==NULL){
        printf("[warning.]Allocate fault !\n");
        return NULL;
    }
    if(initial_value==0)
        return ary;
    else{
        for(i=0; i<ary_size; i++)
            ary[i]=initial_value;
        return ary;
    }
}
//Dynamic 宣告 double array
double* double_Array(int ary_size, double initial_value){
    int i=0;
    double* ary=calloc(ary_size,sizeof(double));
    if(ary==NULL){
        printf("[warning.]Allocate fault !\n");
        return NULL;
    }

    if(initial_value==0)
        return ary;
    else{
        for(i=0; i<ary_size; i++){
            ary[i]=initial_value;
        }
        return ary;
    }
}

void Parameter_Display(char* input_file){
    char buffer[26];    //抓時間用

    printf("*******************[Parameter]*******************\n");
    printf("               %s\n",Get_Time(buffer,1));
    printf("[File]%s\n\n",input_file);
    printf("%d x %d     I=%d    R=%.5f\n",ROW,COLUMN,I_MATRIX,CODE_RATE);
    if(STOP_CRITERION!=0)
    printf("Level-1 block= %d    Level-2 block= %d\n",LEVEL_1_BLK,LEVEL_2_BLK);
    printf("Iteration Max       =%d\n",ITERATION_MAX);
    printf("Stop Criterion mode =%d\n",STOP_CRITERION);
    printf("Erroneous Set Max   =%d\n",ERRONEOUS_SET_MAX);
    printf("SNRdB=%.5f \n",SNRdB);
    printf("*************************************************\n\n");
}

char* Get_Time(char* buffer,int mode){
    // mod1 time display
    // mod2 txt output
    time_t timer;
    struct tm* time_info;

    time(&timer);
    time_info=localtime(&timer);
    if(mode==1){
        strftime(buffer,26,"%Y-%m-%d %H:%M:%S",time_info);
    }else{
        strftime(buffer,40,"backup\\[Backup]%Y-%m-%d-%H-%M.txt",time_info);
    }
    return buffer;
}

void Encoder(int* codeword, int Info_length, int** G_matrix, int G_row, int G_col){
    int i,j;
    int* info_bits=int_Array(Info_length,0);
    int tmp_codeword=0;                 //累加用

    //information bit generate
    for(i=0;i<Info_length;i++){
        info_bits[i]=rand()%2;
    }

    //encoder uG=c
    for(i=0;i<G_row;i++)
    {
        for(j=0;j<G_col;j++)
        {
            if(G_matrix[i][j]==-1)      //到尾端 此 column 沒1了
                break;
            tmp_codeword=tmp_codeword^info_bits[ G_matrix[i][j]-1 ];    //GF(2)
        }
        codeword[i]=tmp_codeword;
        tmp_codeword=0;                 //reset
    }
    free(info_bits);
}
