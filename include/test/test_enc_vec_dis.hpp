#include "seal/seal.h"
#include "seal/util/iterator.h"
#include <numeric>
#include <stdio.h>

using namespace std;
using namespace seal;


void enc_vec_dis_test(){
    cout <<"Task: compute distance between encrypted vectors. "<<endl;
    cout <<"Method: LT(enc(input)*enc(v1||v2||...||vn))"<<endl;
    
    EncryptionParameters parms(scheme_type::bfv);
    size_t poly_modulus_degree = 8192;
    parms.set_poly_modulus_degree(poly_modulus_degree);

    //parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));

    auto bfv_coeff_modulus = CoeffModulus::Create(poly_modulus_degree, { 49,50,60});
    parms.set_coeff_modulus(bfv_coeff_modulus);

    parms.set_plain_modulus(65537);

    //parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20));


    SEALContext context(parms);
    print_parameters(context);
    cout <<endl;

    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    BatchEncoder batch_encoder(context);
    size_t slot_count = batch_encoder.slot_count();
    cout <<"Slot count = "<<slot_count<<endl;

    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    GaloisKeys galois_keys;
    keygen.create_galois_keys(galois_keys);
    cout <<"RLWE sk, pk, Relinearization key and Galois key generated. "<<endl;


    auto &context_data = *context.key_context_data();
    auto coeff_modulus = context_data.parms().coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    uint64_t rlwe_ct_modulus = coeff_modulus[0].value();
    
    //cout <<coeff_modulus[0].value()<<" "<<rlwe_ct_modulus <<endl;
    //cout << coeff_modulus[coeff_modulus_size-1].bit_count() <<endl;

    int vec_len = 512;
    int num_vec = poly_modulus_degree/vec_len;
    int num_batch = 40;
    int pre_store_vec_num = num_vec*num_batch;
    cout <<"vector length = "<<vec_len<<", number of pre stored encrypted vectors = "<<pre_store_vec_num<<endl;
    cout <<"Number of pre stroed vector per batch = "<<num_vec<<endl;

    //generate test vectors
    //input test vector
    vector<int64_t> input_vec(vec_len,0);
    for (int i = 0; i < vec_len; ++i){
        input_vec[i] = 1;
    }

    //pre-stored vectors
    vector<vector<int64_t>> pre_stored_vec(pre_store_vec_num,vector<int64_t>(vec_len,0));
    for (int i = 0; i < pre_store_vec_num; ++i){
        for (int j = 0; j < vec_len; ++j){
        pre_stored_vec[i][j] = 1; 
        }
    }

    cout <<"test vectors generated. "<<endl;

    //prepare encrypted input vector
    vector<uint64_t> pod_matrix(slot_count,0);
    for (int i = 0; i < num_vec; ++i){
        for(int j = 0 ; j < vec_len ; ++j){
            pod_matrix[i*vec_len+j] = input_vec[j];
        }
    }

    Plaintext plain_matrix;
    cout << "input vector encoded and encrypted."<<endl;
    batch_encoder.encode(pod_matrix, plain_matrix);
    Ciphertext encrypted_matrix;
    encryptor.encrypt(plain_matrix, encrypted_matrix);
    cout << "Noise budget in fresh encryption: " << decryptor.invariant_noise_budget(encrypted_matrix) << " bits"<< endl;
    cout <<"number of primes in coeff_modulus of ct: "<<encrypted_matrix.coeff_modulus_size()<<endl;

    //evaluator.transform_to_ntt_inplace(encrypted_matrix);

    //cout <<"Encrypted input vector transformed to ntt space. "<<endl;

    //prepare encrypted pre sotred vector
    vector<Ciphertext> encrypt_pre_store_vector(num_batch);

    for (int i = 0; i < num_batch; ++i){
        vector<uint64_t> pod_matrix_2(slot_count,0);
        for (int k = 0; k < num_vec; ++k){
            for(int j = 0 ; j < vec_len ; ++j){
                pod_matrix_2[k*vec_len+j] = pre_stored_vec[i*num_vec+k][j];
            }
        }

        Plaintext plain_matrix_2;
        //cout << "Encode and encrypt.";
        batch_encoder.encode(pod_matrix_2, plain_matrix_2);
        Ciphertext encrypted_matrix_2;
        encryptor.encrypt(plain_matrix_2, encrypted_matrix_2);
        encrypt_pre_store_vector[i] = encrypted_matrix_2;

    }

    cout << "pre-sotred vectors encoded and encrypted. Number of ciphertext = "<<encrypt_pre_store_vector.size()<<endl;

    //prepare LT matrix
    vector<vector<uint64_t>> M(poly_modulus_degree, vector<uint64_t>(poly_modulus_degree,0));
    for (int i = 0; i < num_vec; ++i){
        for(int j = 0 ; j < vec_len ; ++j){
            for(int k = 0 ; k < vec_len ; ++k){
                M[i*vec_len+j][i*vec_len+k] = 1;
            }
        }
    }

    vector<Plaintext> out = ecd_M(M,poly_modulus_degree,context);

    cout <<"LT matrix generated and encoded. "<<endl;


    struct timeval tstart, tend;

    cout<<"----------------------TEST START-----------------------"<<endl;
    gettimeofday(&tstart,NULL);

    vector<Ciphertext> eval_result(num_batch);

    #pragma omp parallel for

    for (int i = 0; i < num_batch; ++i){
        //cout <<"Batch No. "<<i+1<<endl;

        Ciphertext mult;
        evaluator.multiply(encrypted_matrix, encrypt_pre_store_vector[i], mult);
        evaluator.relinearize_inplace(mult, relin_keys);
        //evaluator.mod_switch_to_next_inplace(mult);
        //cout <<"number of primes in coeff_modulus of ct: "<<mult.coeff_modulus_size()<<endl;
        

        Ciphertext encrypted_result = LT_ecd_M(mult, out, poly_modulus_degree, galois_keys, context,secret_key);
       //cout << "Noise budget in LT result: " << decryptor.invariant_noise_budget(encrypted_result) << " bits"<< endl;
        eval_result[i] = encrypted_result;

        
    }

    gettimeofday(&tend,NULL);

    cout <<"-------------------------TEST END--------------------"<<endl;

    double  time_ek = tend.tv_sec-tstart.tv_sec+(tend.tv_usec-tstart.tv_usec)/1000000.0;
    cout <<"Total time for evaluations: "<<time_ek<<"s. "<<endl;
    cout <<"Average time for per evaluation: "<<time_ek/(double)pre_store_vec_num<<"s. "<<endl;

    for (int i = 0; i < num_batch; ++i){
        Plaintext plain_result;
        decryptor.decrypt(eval_result[i], plain_result);
        vector<uint64_t> dcd_matrix(slot_count,0);
        batch_encoder.decode(plain_result, dcd_matrix);
        cout<<"Test result: ";
        for (int i = 0; i < poly_modulus_degree; ++i){
            if(i % vec_len == 0){
                cout <<dcd_matrix[i]<<" ";
            }
            
        }
        cout <<endl;
    }

    cout <<endl;





    
    

}