#include "seal/seal.h"
#include "seal/util/iterator.h"
#include <numeric>
#include <stdio.h>

using namespace std;
using namespace seal;

inline void print_parameters(const seal::SEALContext &context)
{
    auto &context_data = *context.key_context_data();

    /*
    Which scheme are we using?
    */
    std::string scheme_name;
    switch (context_data.parms().scheme())
    {
    case seal::scheme_type::bfv:
        scheme_name = "BFV";
        break;
    case seal::scheme_type::ckks:
        scheme_name = "CKKS";
        break;
    case seal::scheme_type::bgv:
        scheme_name = "BGV";
        break;
    default:
        throw std::invalid_argument("unsupported scheme");
    }
    std::cout << "/" << std::endl;
    std::cout << "| Encryption parameters :" << std::endl;
    std::cout << "|   scheme: " << scheme_name << std::endl;
    std::cout << "|   poly_modulus_degree: " << context_data.parms().poly_modulus_degree() << std::endl;

    /*
    Print the size of the true (product) coefficient modulus.
    */
    std::cout << "|   coeff_modulus size: ";
    std::cout << context_data.total_coeff_modulus_bit_count() << " (";
    auto coeff_modulus = context_data.parms().coeff_modulus();
    std::size_t coeff_modulus_size = coeff_modulus.size();
    for (std::size_t i = 0; i < coeff_modulus_size - 1; i++)
    {
        std::cout << coeff_modulus[i].bit_count() << " + ";
    }
    std::cout << coeff_modulus.back().bit_count();
    std::cout << ") bits" << std::endl;

    /*
    For the BFV scheme print the plain_modulus parameter.
    */
    if (context_data.parms().scheme() == seal::scheme_type::bfv)
    {
        std::cout << "|   plain_modulus: " << context_data.parms().plain_modulus().value() << std::endl;
    }

    std::cout << "\\" << std::endl;
}


void enc_vec_dis_batch_test(){
    cout <<"Task: compute distance between encrypted vectors. "<<endl;
    cout <<"Method: SUM(Enc(input_i)*Enc(v1i,v2i,...,vni))"<<endl;
    
    EncryptionParameters parms(scheme_type::bfv);
    size_t poly_modulus_degree = 8192;
    parms.set_poly_modulus_degree(poly_modulus_degree);

    parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));

    //auto bfv_coeff_modulus = CoeffModulus::Create(poly_modulus_degree, { 44,44,44});
    //parms.set_coeff_modulus(bfv_coeff_modulus);

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

    /*
    GaloisKeys galois_keys;
    keygen.create_galois_keys(galois_keys);
    */
    cout <<"RLWE sk, pk, Relinearization key generated. "<<endl;


    auto &context_data = *context.key_context_data();
    auto coeff_modulus = context_data.parms().coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    uint64_t rlwe_ct_modulus = coeff_modulus[0].value();
    
    //cout <<coeff_modulus[0].value()<<" "<<rlwe_ct_modulus <<endl;
    //cout << coeff_modulus[coeff_modulus_size-1].bit_count() <<endl;

    int vec_len = 768;
    int num_vec = poly_modulus_degree;
    int num_batch = 1;
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
    vector<Ciphertext> enc_input(vec_len);

    for (int i = 0; i < vec_len; ++i){
        vector<uint64_t> pod_matrix(slot_count,input_vec[i]);
        Plaintext plain_matrix;
        batch_encoder.encode(pod_matrix, plain_matrix);
        Ciphertext encrypted_matrix;
        encryptor.encrypt(plain_matrix, encrypted_matrix);
        enc_input[i] = encrypted_matrix;
        //evaluator.transform_to_ntt_inplace(enc_input[i]);

    //cout <<"Encrypted input vector transformed to ntt space. "<<endl;
    }

    cout <<"Input vector is encoded and encrypted. "<<endl;
    cout << "Noise budget in fresh encryption: " << decryptor.invariant_noise_budget(enc_input[0]) << " bits"<< endl;
    cout <<"number of primes in coeff_modulus of ct: "<<enc_input[0].coeff_modulus_size()<<endl;

    //prepare encrypted pre sotred vector
    vector<Ciphertext> encrypt_pre_store_vector(vec_len);

    for (int i = 0; i < vec_len; ++i){
        vector<uint64_t> pod_matrix_2(slot_count,0);
        for(int j = 0 ; j < slot_count ; ++j){
                pod_matrix_2[j] = pre_stored_vec[j][i];
        }

        Plaintext plain_matrix_2;
        //cout << "Encode and encrypt.";
        batch_encoder.encode(pod_matrix_2, plain_matrix_2);
        Ciphertext encrypted_matrix_2;
        encryptor.encrypt(plain_matrix_2, encrypted_matrix_2);
        encrypt_pre_store_vector[i] = encrypted_matrix_2;
        //evaluator.transform_to_ntt_inplace(encrypt_pre_store_vector[i]);

    }

    cout << "pre-sotred vectors encoded and encrypted. Number of ciphertext = "<<encrypt_pre_store_vector.size()<<endl;


    struct timeval tstart, tend;

    cout<<"----------------------TEST START-----------------------"<<endl;
    gettimeofday(&tstart,NULL);

    vector<Ciphertext> eval_result(vec_len);

    #pragma omp parallel for

    for (int i = 0; i < vec_len; ++i){
        //cout <<"Batch No. "<<i+1<<endl;

        Ciphertext mult;
        evaluator.multiply(enc_input[i], encrypt_pre_store_vector[i], mult);
        evaluator.relinearize_inplace(mult, relin_keys);
        //evaluator.mod_switch_to_next_inplace(mult);
        //cout <<"number of primes in coeff_modulus of ct: "<<mult.coeff_modulus_size()<<endl;
        eval_result[i] = mult;
    }

    Ciphertext add_result = eval_result[0];

    for (int i = 1; i < vec_len; ++i){
        evaluator.add_inplace(add_result,eval_result[i]);
    }

    gettimeofday(&tend,NULL);

    cout <<"-------------------------TEST END--------------------"<<endl;

    double  time_ek = tend.tv_sec-tstart.tv_sec+(tend.tv_usec-tstart.tv_usec)/1000000.0;
    cout <<"Total time for evaluations: "<<time_ek<<"s. "<<endl;
    cout <<"Average time for per evaluation: "<<time_ek/(double)pre_store_vec_num<<"s. "<<endl;

    Plaintext plain_result;
    decryptor.decrypt(add_result, plain_result);
    vector<uint64_t> dcd_matrix(slot_count,0);
    batch_encoder.decode(plain_result, dcd_matrix);
    cout<<"Test result of first 10 ct: ";
    for (int i = 0; i < 10; ++i){
        cout <<dcd_matrix[i]<<" ";
    }    
    cout <<endl;
    cout <<endl;


}