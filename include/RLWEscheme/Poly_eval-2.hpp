using namespace std;
using namespace seal;

inline std::string uint64_to_hex_string(std::uint64_t value)
{
    return seal::util::uint_to_hex_string(&value, std::size_t(1));
}

//compute RLWE(ecd(f(x))) from RLWE(ecd(x)) and f(x)
extern Ciphertext poly_eval(const Ciphertext & RLWE_ct, 
const RelinKeys & RelK, const SEALContext& seal_context, const vector<uint64_t> & f,
const SecretKey & sk){

  Evaluator evaluator(seal_context);

  //for test
  Decryptor decryptor(seal_context, sk);
  BatchEncoder batch_encoder(seal_context);

  int f_degree = f.size()-1;
  Ciphertext out(RLWE_ct);

  Plaintext pt_f_n(uint64_to_hex_string(f[f_degree]));
  evaluator.multiply_plain_inplace(out,pt_f_n);

  Plaintext pt_f_n_minus_1(uint64_to_hex_string(f[f_degree-1]));
  evaluator.add_plain_inplace(out,pt_f_n_minus_1);

/*
  cout <<"i = "<<f_degree-1<<endl;
    Plaintext test_result;
    decryptor.decrypt(out,test_result);
    vector<uint64_t> dcd_test;
    batch_encoder.decode(test_result,dcd_test);
    for(int j = 0 ; j < 10 ; ++j){
      cout <<dcd_test[j]<<" ";
    }
    cout <<"...... ";
    for (int j = dcd_test.size()-10; j < dcd_test.size(); ++j){
      cout <<dcd_test[j]<<" ";
    }
    cout <<endl;

  */

  for(int i = f_degree-2 ; i >= 0 ; --i){
    evaluator.multiply_inplace(out,RLWE_ct);
    evaluator.relinearize_inplace(out,RelK);
    Plaintext pt_f_i(uint64_to_hex_string(f[i]));
    evaluator.add_plain_inplace(out,pt_f_i);

/*
     //test correctness 
    cout <<"i = "<<i<<endl;
    Plaintext test_result;
    decryptor.decrypt(out,test_result);
    vector<uint64_t> dcd_test;
    batch_encoder.decode(test_result,dcd_test);
    for(int j = 0 ; j < 10 ; ++j){
      cout <<dcd_test[j]<<" ";
    }
    cout <<"...... ";
    for (int j = dcd_test.size()-10; j < dcd_test.size(); ++j){
      cout <<dcd_test[j]<<" ";
    }
    cout <<endl;
    cout << "Noise budget in poly_eval result: " << decryptor.invariant_noise_budget(out) << " bits"<< endl;
    */

  }
  return out;
}

inline void Parallel_calUptoDegreeK(vector<Ciphertext>& output, const Ciphertext& input, const int DegreeK, const RelinKeys &relin_keys,
                           const SEALContext& context, const bool skip_odd=false){

    int num_cal = 0;

    vector<int> calculated(DegreeK, 0); // 0 is for degree 1, 1 is for degree 2
    calculated[0] = 1;
    Evaluator evaluator(context);
    output[0] = input;
    num_cal ++;

    vector<vector<int>> combined(DegreeK,vector<int>(DegreeK,0));

/*
    //skip odd degree if skip_odd = true
    for(int i = DegreeK; i > 0; --i){
        if (skip_odd && i % 2 == 1) { // 0 is for degree 1, 1 is for degree 2, skip all 2k+1 degree
            calculated[i-1] = 1;
            output[i-1] = input;
        }
    }
*/

    //step 1, compute x^2, x^4, x^8,....in sequence
    int index = 1;
    int num_cal_1 = 1;
  //  cout <<"calculated in step 1: "<<endl;
    while(index < DegreeK){
        int new_index = index * 2;
        Ciphertext temp = output[index-1];
        evaluator.square_inplace(temp);
        evaluator.relinearize_inplace(temp, relin_keys);
        //cout <<"square end. "<<endl;
        output[new_index-1] = temp;
        calculated[new_index-1] = 1;
        num_cal_1++;
       // cout <<new_index<<" ";
        index = new_index;
    }
  //  cout <<endl;

    //step 2, compute x^{1,2,4,...,}*x^{1,2,4,...}
   // cout <<"calculated in step 2: "<<endl;

    #pragma omp parallel for

    for(int ii = 0 ; ii < num_cal_1 ; ++ii){
        int i = 1 <<ii;
        //cout <<"i = "<<i<<endl;
        for(int j = i*2 ; j < DegreeK ; j *= 2){
        //    cout <<"j = "<<j<<", i+j = "<<i+j<<endl;
            if(i + j > DegreeK){
                break;
            }
            combined[i-1][j-1] = 1;
            int new_index = i + j;
            if (i == j){
                continue;
            }
            if (calculated[new_index-1] == 1){
               continue;
            }
            else{
                Ciphertext temp = output[i-1];
                Ciphertext resct = output[j-1];
                evaluator.multiply_inplace(temp, resct);
                evaluator.relinearize_inplace(temp, relin_keys);
                output[new_index-1] = temp;
                calculated[new_index-1] = 1;
            //    cout <<"calculated = "<<new_index <<endl;
            }
        }
    }
  //  cout <<endl;

    if(skip_odd == 1){
        cout <<"to be finished. "<<endl;
    }

    else{

        #pragma omp parallel for

        for (int i = 1; i < DegreeK; ++i){

            if(calculated[i-1] == 0){
                continue;
            }
            for (int j = i+1; j < DegreeK; ++j){

                if(i + j > DegreeK){
                break;
                }


                if(calculated[j-1] == 0){
                    continue;
                }

                if(combined[i-1][j-1] == 1){
                    continue;
                }

                int new_index = i + j;
                combined[i-1][j-1] = 1;

                if (calculated[new_index-1] == 1){
                    continue;
                }
                else{
                    Ciphertext temp = output[i-1];
                    Ciphertext resct = output[j-1];
                    evaluator.multiply_inplace(temp, resct);
                    evaluator.relinearize_inplace(temp, relin_keys);
                    output[new_index-1] = temp;
                    calculated[new_index-1] = 1;
                }

            }
        }
    }

/*
    //check how many ct are calculated
    
    for (int i = 0; i < DegreeK; ++i){
        if(calculated[i] == 1){
            cout <<i+1 <<" ";
            num_cal ++;
        }
    }
    cout <<endl;
    cout <<"number of calculated ct = "<<num_cal<<endl;
*/
    int out_size = output.size();

    for(size_t i = 0; i < out_size-1; i++){
        evaluator.mod_switch_to_inplace(output[i], output[out_size-1].parms_id()); // match modulus
    }
    return;


}


//copy from asiacrypt23

inline void calUptoDegreeK(vector<Ciphertext>& output, const Ciphertext& input, const int DegreeK, const RelinKeys &relin_keys,
                           const SEALContext& context, const bool skip_odd=false) {
    vector<int> calculated(DegreeK, 0);
    Evaluator evaluator(context);
    output[0] = input;
    calculated[0] = 1; // degree 1, x
    Ciphertext res, base;
    vector<int> numMod(DegreeK, 0);

    int mod_switch_num = 0;

    for(int i = DegreeK; i > 0; --i){
        if (skip_odd && i % 2 == 1) { // 0 is for degree 1, 1 is for degree 2, skip all 2k+1 degree
            calculated[i-1] = 1;
            //cout <<i <<" ";
            output[i-1] = input;
        }
        else if(calculated[i-1] == 0){
            auto toCalculate = i;
            int resdeg = 0;
            int basedeg = 1;
            base = input;
            while(toCalculate > 0){
                if(toCalculate & 1){
                    toCalculate -= 1;
                    resdeg += basedeg;
                    if(calculated[resdeg-1] != 0){
                        res = output[resdeg - 1];
                    } 
                    else {
                        if(resdeg == basedeg){
                            res = base; // should've never be used as base should have made calculated[basedeg-1]
                        }
                        else {
                            numMod[resdeg-1] = numMod[basedeg-1];

                            evaluator.mod_switch_to_inplace(res, base.parms_id()); // match modulus
                            evaluator.multiply_inplace(res, base);
                            evaluator.relinearize_inplace(res, relin_keys);
                            while(numMod[resdeg-1] < (ceil(log2(resdeg))/2)){
                               // cout <<resdeg<<" ";
                                evaluator.mod_switch_to_next_inplace(res);
                                mod_switch_num++;
                                numMod[resdeg-1]+=1;
                            }
                        }
                        output[resdeg-1] = res;
                        calculated[resdeg-1] += 1;
                        cout <<resdeg<<" ";
                    }
                } 
                else {
                    toCalculate /= 2;
                    basedeg *= 2;
                    if(calculated[basedeg-1] != 0){
                        base = output[basedeg - 1];
                    }
                    else {
                        numMod[basedeg-1] = numMod[basedeg/2-1];
                        evaluator.square_inplace(base);
                        evaluator.relinearize_inplace(base, relin_keys);
                        while(numMod[basedeg-1] < (ceil(log2(basedeg))/2)){
                           // cout <<basedeg<<" ";
                            evaluator.mod_switch_to_next_inplace(base);
                            mod_switch_num++;
                            numMod[basedeg-1]+=1;
                        }
                        output[basedeg-1] = base;
                        calculated[basedeg-1] += 1;
                        cout <<basedeg<<" ";
                    }
                }
            }
        }
    }

    cout<<mod_switch_num<<endl;
    for(size_t i = 0; i < output.size()-1; i++){
        evaluator.mod_switch_to_inplace(output[i], output[output.size()-1].parms_id()); // match modulus
    }
    return;
}


void Bootstrap_RangeCheck_PatersonStockmeyer(Ciphertext& ciphertext, const Ciphertext& input, const vector<uint64_t>& rangeCheckIndices,
                                            const size_t& degree, const RelinKeys &relin_keys, const SEALContext& context, const SecretKey& bfv_secret_key, 
                                             const int f_zero = 0, const bool skip_first_odd = false, const int firstDegree = 256, const int secondDegree = 256) {
    MemoryPoolHandle my_pool_larger = MemoryPoolHandle::New(true);
    auto old_prof_larger = MemoryManager::SwitchProfile(std::make_unique<MMProfFixed>(std::move(my_pool_larger)));

    Decryptor decryptor(context, bfv_secret_key);
    Evaluator evaluator(context);
    BatchEncoder batch_encoder(context);

    struct timeval tstart1, tend1;

    gettimeofday(&tstart1,NULL);

    vector<Ciphertext> kCTs(firstDegree);
    Parallel_calUptoDegreeK(kCTs, input, firstDegree, relin_keys, context,skip_first_odd);


    vector<Ciphertext> kToMCTs(secondDegree);
    Parallel_calUptoDegreeK(kToMCTs, kCTs[kCTs.size()-1], secondDegree, relin_keys, context, skip_first_odd);

    for (int j = 0; j < (int) kCTs.size(); j++) {
      evaluator.mod_switch_to_inplace(kCTs[j], kToMCTs[kToMCTs.size()-1].parms_id());
    }

    for (int j = 0; j < (int) kToMCTs.size(); j++) {
      evaluator.mod_switch_to_next_inplace(kToMCTs[j]);
    }
    cout << "Noise for last: " << decryptor.invariant_noise_budget(kToMCTs[kToMCTs.size()-1]) << " bits\n";

    gettimeofday(&tend1,NULL);
    double new_x_k_time = tend1.tv_sec-tstart1.tv_sec+(tend1.tv_usec-tstart1.tv_usec)/1000000.0;
    cout <<"Time for Parallel_calUptoDegreeK = "<<new_x_k_time<<endl;

/*
    size_t slot_count = batch_encoder.slot_count();
    Plaintext plain_result;
    decryptor.decrypt(kToMCTs[secondDegree-2], plain_result);
    vector<uint64_t> dcd_matrix(slot_count,0);
    batch_encoder.decode(plain_result, dcd_matrix);

    for(int j = 0 ; j < 10 ; ++j){
      cout <<dcd_matrix[j]<<" ";
    }
    cout <<"...... ";
    for (int j = dcd_matrix.size()-10; j < dcd_matrix.size(); ++j){
      cout <<dcd_matrix[j]<<" ";
    }
    cout <<endl;
*/

/*
    gettimeofday(&tstart1,NULL);

    calUptoDegreeK(kCTs, input, firstDegree, relin_keys, context, skip_first_odd);

    //vector<Ciphertext> kToMCTs(secondDegree);
    calUptoDegreeK(kToMCTs, kCTs[kCTs.size()-1], secondDegree, relin_keys, context);

    for (int j = 0; j < (int) kCTs.size(); j++) {
      evaluator.mod_switch_to_inplace(kCTs[j], kToMCTs[kToMCTs.size()-1].parms_id());
    }

    for (int j = 0; j < (int) kToMCTs.size(); j++) {
      evaluator.mod_switch_to_next_inplace(kToMCTs[j]);
    }
    cout << "Noise for last: " << decryptor.invariant_noise_budget(kToMCTs[kToMCTs.size()-1]) << " bits\n";

    gettimeofday(&tend1,NULL);
    double x_k_time = tend1.tv_sec-tstart1.tv_sec+(tend1.tv_usec-tstart1.tv_usec)/1000000.0;
    cout <<"Time for calUptoDegreeK = "<<x_k_time<<endl;
*/
    

    gettimeofday(&tstart1,NULL);

    Ciphertext temp_relin;
    vector<Ciphertext> levelSum(secondDegree);
    // int third = 0;
    
    Plaintext plainInd;
    plainInd.resize(degree);
    plainInd.parms_id() = parms_id_zero;
    for (int i = 0; i < (int) degree; i++) {
        plainInd.data()[i] = 0;
    }
   // cout <<"init plainInd"<<endl;

    #pragma omp parallel for

    for(int i = 0; i < secondDegree; i++) {
      //  Ciphertext levelSum;
        bool flag = false;
        for(int j = 0; j < firstDegree; j++) {
         // cout <<i*firstDegree+j<<" f: "<<rangeCheckIndices[i*firstDegree+j]<<endl;
            if(rangeCheckIndices[i*firstDegree+j] != 0) {
                // vector<uint64_t> intInd(degree, rangeCheckIndices[i*firstDegree+j]);
                plainInd.data()[0] = rangeCheckIndices[i*firstDegree+j];
                // batch_encoder.encode(intInd, plainInd);
                if (!flag) {
                    evaluator.multiply_plain(kCTs[j], plainInd, levelSum[i]);
                    flag = true;
                }
                else {
                    Ciphertext tmp;
                    evaluator.multiply_plain(kCTs[j], plainInd, tmp);
                    evaluator.add_inplace(levelSum[i], tmp);
                }
            }
        }
        evaluator.mod_switch_to_inplace(levelSum[i], kToMCTs[i].parms_id()); // mod down the plaintext multiplication noise
        // time_start = chrono::high_resolution_clock::now();
      }

      for (int i = 0; i < secondDegree; ++i){
        if(i == 0) {
            ciphertext = levelSum[0];
        }
        else if (i == 1) { // initialize for temp_relin, which is of ct size = 3
            evaluator.multiply(levelSum[i], kToMCTs[i - 1], temp_relin);
        }
        else {
            evaluator.multiply_inplace(levelSum[i], kToMCTs[i - 1]);
            evaluator.add_inplace(temp_relin, levelSum[i]);
        }
        // time_end = chrono::high_resolution_clock::now();
        // third += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    }
    gettimeofday(&tend1,NULL);
    double eval_time = tend1.tv_sec-tstart1.tv_sec+(tend1.tv_usec-tstart1.tv_usec)/1000000.0;
    cout <<"Time for coeff*x^k = "<<eval_time<<endl;
  
    /*
    for(int i = 0; i < firstDegree; i++){
        kCTs[i].release();
    }
    for(int i = 0; i < secondDegree; i++){
        kToMCTs[i].release();
    }
*/

    evaluator.relinearize_inplace(temp_relin, relin_keys);
    evaluator.add_inplace(ciphertext, temp_relin);
  //  temp_relin.release();

    // vector<uint64_t> intInd(degree, f_zero); 
    // Plaintext plainInd;
    //plainInd.data()[0] = f_zero;
    // batch_encoder.encode(intInd, plainInd);
   // evaluator.negate_inplace(ciphertext);
   // evaluator.add_plain_inplace(ciphertext, plainInd);

    cout << "Noise after function: " << decryptor.invariant_noise_budget(ciphertext) << " bits\n";
    
    MemoryManager::SwitchProfile(std::move(old_prof_larger));
}
