# -*- coding: utf-8 -*-
# %tensorflow_version 2.x
import os
import numpy as np
import tensorflow as tf

print("Tensorflow version " + tf.__version__)

import pandas as pd
from keras.models import load_model
import time
import pickle as pk

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 50)
pd.set_option('display.width', 1000)


def _error(actual: np.ndarray, predicted: np.ndarray):
    """ Simple error """
    return actual - predicted


def bias(actual: np.ndarray, predicted: np.ndarray):
    """ BIAS """
    return np.mean(_error(actual, predicted))


def mse(actual: np.ndarray, predicted: np.ndarray):
    """ Mean Squared Error """
    return np.mean(np.square(_error(actual, predicted)))


def rmse(actual: np.ndarray, predicted: np.ndarray):
    """ Root Mean Squared Error """
    return np.sqrt(mse(actual, predicted))


def mae(actual: np.ndarray, predicted: np.ndarray):
    """ Mean Absolute Error """
    return np.mean(np.abs(_error(actual, predicted)))


def _percentage_error(actual: np.ndarray, predicted: np.ndarray):
    """
    Percentage error
    Note: result is NOT multiplied by 100
    """
    return _error(actual, predicted) / (actual + EPSILON)


def mpe(actual: np.ndarray, predicted: np.ndarray):
    """ Mean Percentage Error """
    return np.mean(_percentage_error(actual, predicted))


def mape(actual: np.ndarray, predicted: np.ndarray):
    """
    Mean Absolute Percentage Error
    Properties:
        + Easy to interpret
        + Scale independent
        - Biased, not symmetric
        - Undefined when actual[t] == 0
    Note: result is NOT multiplied by 100
    """
    return np.mean(np.abs(_percentage_error(actual, predicted)))


def rmspe(actual: np.ndarray, predicted: np.ndarray):
    """
    Root Mean Squared Percentage Error
    Note: result is NOT multiplied by 100
    """
    return np.sqrt(np.mean(np.square(_percentage_error(actual, predicted))))


def preprocess_features(input_df):
    """Prepares input features input_df

  Args:
    input data frame: A Pandas DataFrame expected to contain data from input data set.
  Returns:
    A DataFrame that contains the features to be used for the model, including
    synthetic features.
  """
    # All input variables
    # "k", "si", "si_kmax_+_1", "sl", "Tc", "qv", "qc", "qr", "qi", "qs", "qg", "ni", "ns", "nr", "NG", "NC", "tke", "kzh", "gps", "omega"
    input_features = ["k", "sl", "Tc", "qv", "qc", "qr", "qi", "qs", "qg", "ni", "ns", "nr", "NG", "NC", "tke", "omega",
                      "LSRAIN", "LSSNOW"]
    selected_features = input_df[input_features]

    # processed_features = selected_features.copy()
    return selected_features, input_features


def preprocess_targets(output_df):
    """Prepares target features (i.e., labels) from output_df

  Args:
    output_df: A Pandas DataFrame expected to contain data from output data set.
  Returns:
    A DataFrame that contains the target features.
  """
    # All output variables
    # "k", "Tc", "qv", "qc", "qr", "qi", "qs", "qg", "ni", "ns", "nr", "NG", "NC", "EFFCS", "EFFIS", "LSRAIN", "LSSNOW"
    output_targets = output_df[
        ["k", "Tc", "qv", "qc", "qr", "qi", "qs", "qg", "ni", "ns", "nr", "NG", "NC", "LSRAIN", "LSSNOW"]]
    return output_targets


def linear_scale(serie_or_np_arr, min_val, max_val):
    # https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.MinMaxScaler.html
    #   X_std = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))
    # X_scaled = X_std * (max - min) + min
    min, max = 0, 1
    y = lambda x: ((x - min_val) / (max_val - min_val)) * (max - min) + min
    return y(serie_or_np_arr)
    # Linear normalization (serie or numpy 1D array)
    # -1 e 1


def delinear_scale(serie_or_np_arr, min_val, max_val):
    # https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.MinMaxScaler.html
    #   X_std = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))
    # X_scaled = X_std * (max - min) + min
    min, max = 0, 1
    y = lambda x: min_val + (max_val - min_val) * (x - min) / (max - min)
    return y(serie_or_np_arr)


def normalize_linear_scale(examples_dataframe, minmax_dict=None, scaler="normalize", minmax_border_perc=0.0):
    """Returns a version of the input `DataFrame` that has all its features normalized linearly."""

    # Convert pandas data into a dict of np arrays.
    processed_features = examples_dataframe.copy()

    for key, value in dict(examples_dataframe).items():
        if key != 'k':
            if minmax_dict is None:
                min_val = value.min() - value.min() * minmax_border_perc
                max_val = value.max() + value.max() * minmax_border_perc
            else:
                min_val = minmax_dict[key][0]
                max_val = minmax_dict[key][1]
            if scaler == "normalize":
                processed_features[key] = linear_scale(value, min_val, max_val)
            # else:
            #   processed_features[key] = linear_scale_htang(value, min_val, max_val)

    return processed_features


def denormalize_linear_scale(examples_dataframe, minmax_dict=None, scaler="normalize", minmax_border_perc=0.0,
                             is_col_k=False):
    """Returns a version of the input `DataFrame` that has all its features normalized linearly."""

    # Convert pandas data into a dict of np arrays.
    processed_features = examples_dataframe.copy()

    for key, value in dict(examples_dataframe).items():
        if key != 'k':
            if minmax_dict is None:
                min_val = value.min() - value.min() * minmax_border_perc
                max_val = value.max() + value.max() * minmax_border_perc
            else:
                if is_col_k:
                    col_name_k = key.split("_k")
                    key_minmax = col_name_k[0]
                else:
                    key_minmax = key
                min_val = minmax_dict[key_minmax][0]
                max_val = minmax_dict[key_minmax][1]
            if scaler == "normalize":
                processed_features[key] = delinear_scale(value, min_val, max_val)
            # else:
            #   processed_features[key] = delinear_scale_htang(value, min_val, max_val)

    return processed_features


def get_dic_levels_excluded(use_levs):
    if use_levs == 'Lev':
        if k_max == 18:  # 2014-2015 dt 360 ensemble 5x
            dic_var_levels_exclude_examples = { \
                'ni': list(range(1, 7)) + list(range(16, 19)),
                'ns': [17, 18],
                'nr': range(12, 19),
                'NG': [17, 18],
                'NC': range(12, 19),
                'tke': [9, 13, 14, 15, 16, 18],
                'omega': [17, 18]
            }

            dic_var_levels_exclude_targets = { \
                'qc': range(12, 19),
                'qr': range(12, 19),
                'qi': list(range(1, 7)) + list(range(16, 19)),
                'qs': [17, 18],
                'qg': [17, 18],
                'ni': list(range(1, 7)) + list(range(16, 19)),
                'ns': [17, 18],
                'nr': range(12, 19),
                'NG': [17, 18],
                'NC': range(12, 19)
            }
        elif k_max == 28:  # 2014-2015
            dic_var_levels_exclude_examples = { \
                'qc': range(16, 29),
                'qr': range(16, 29),
                'qi': list(range(1, 11)) + list(range(20, 29)),
                'qs': range(23, 29),
                'qg': range(23, 29),
                'ni': list(range(1, 11)) + list(range(20, 29)),
                'ns': list(range(1, 3)) + list(range(23, 29)),
                'nr': range(16, 29),
                'NG': range(23, 29),
                'NC': range(16, 29),
                'omega': range(23, 28),
            }

            dic_var_levels_exclude_targets = { \
                'qc': range(16, 29),
                'qr': range(16, 29),
                'qi': list(range(1, 11)) + list(range(20, 29)),
                'qs': range(23, 29),
                'qg': range(23, 29),
                'ni': list(range(1, 11)) + list(range(20, 29)),
                'ns': list(range(1, 3)) + list(range(23, 29)),
                'nr': range(16, 29),
                'NG': range(23, 29),
                'NC': range(16, 29)
            }

    else:
        dic_var_levels_exclude_examples, dic_var_levels_exclude_targets = None, None

    return dic_var_levels_exclude_examples, dic_var_levels_exclude_targets


def get_dic_default_values():
    # "Tc", "qv", "qc", "qr", "qi", "qs", "qg", "ni", "ns", "nr", "NG", "NC", "EFFCS", "EFFIS", "LSRAIN", "LSSNOW"]

    # 28 lev - confirmado anteriormente ???
    # 'EFFCS': 0.573833,
    # 'EFFIS': 0.119181,

    # 18 lev - confirmando com valor 25
    # EFFCS = [-5.8909090909, 42.3454545454]
    # EFFIS = [-10.400000000000002, 153.4]
    # print(linear_scale(25, EFFCS[0], EFFCS[1]))
    # print(linear_scale(25, EFFIS[0], EFFIS[1]))
    # >> 0.6404070863180329,
    # >> 0.21611721611721613,

    # if k_max == 18:
    dic_var_default_values = { \
            'Tc': 0.0,
            'qv': 0.0,
            'qc': 0.0,
            'qr': 0.0,
            'qi': 0.0,
            'qs': 0.0,
            'qg': 0.0,
            'ni': 0.0,
            'ns': 0.0,
            'nr': 0.0,
            'NG': 0.0,
            'NC': 0.0,
            'tke': 0.0,
            'kzh': 0.0,
            'omega': 0.0,
            'EFFCS': 25.0,
            'EFFIS': 25.0,
            'LSRAIN': 0.0,
            'LSSNOW': 0.0
    }

    return dic_var_default_values



# versão atual que grava primeiro por variavel depois por nivel, e remove niveis de LSRAIN e LSSNOW
def get_df_col_k(df_orig, k_inicial, k_final, dic_var_levs_exclude=None):
    print(df_orig.columns)
    dic_var_levs_exclude_ok = { \
        'LSRAIN': range(2, k_final+1),
        'LSSNOW': range(2, k_final+1)
    }
    if dic_var_levs_exclude is not None:
        dic_var_levs_exclude_ok.update(dic_var_levs_exclude)

    df_col_k = pd.DataFrame()
    for col in df_orig.columns.tolist():
        if col == 'k':
            continue
        for k in range(k_inicial, k_final + 1):
            df = df_orig.copy()
            df_k = df[df['k'] == k]
            if col in dic_var_levs_exclude_ok.keys() and k in dic_var_levs_exclude_ok[col]:
                continue
            df_col_k['{}_k{}'.format(col, k)] = df_k[col].to_numpy()
    return df_col_k


# versão atual que grava primeiro por variavel depois por nivel, e remove niveis de LSRAIN e LSSNOW
# para 1 linha por nível, usar k_final = k_inicial
#
def get_dffull_from_npcol_k(np_col_k, all_cols, k_inicial, k_final, dic_var_levs_exclude=None,
                            dic_defatult_values=None):
    dic_var_levs_exclude_ok = { \
        'LSRAIN': range(2, k_final+1),
        'LSSNOW': range(2, k_final+1)
    }
    if dic_var_levs_exclude is not None:
        dic_var_levs_exclude_ok.update(dic_var_levs_exclude)

    pos = 0
    df_all_cols = pd.DataFrame()
    for k in range(k_inicial, k_final + 1):
        df_all_cols_tmp = pd.DataFrame()
        # col_k = '{}_k{}'.format(col, k)
        for col in all_cols:
            if col == 'k':
                df_all_cols_tmp[col] = pd.Series(k)
                continue

            if col in dic_var_levs_exclude_ok.keys() and k in dic_var_levs_exclude_ok[col]:
                df_all_cols_tmp[col] = pd.Series(dic_defatult_values[col])
            else:
                df_all_cols_tmp[col] = pd.Series(np_col_k[pos])
                pos += 1
        df_all_cols = pd.concat([df_all_cols, df_all_cols_tmp])

    return df_all_cols


# versão atual que grava primeiro por variavel depois por nivel, e remove niveis de LSRAIN e LSSNOW
def set_default_values(df_all_cols, k_inicial, k_final, dic_var_levs_exclude=None, dic_defatult_values=None):
    dic_var_levs_exclude_ok = { \
        'LSRAIN': range(2, k_max),
        'LSSNOW': range(2, k_max)
    }
    if dic_var_levs_exclude is not None:
        dic_var_levs_exclude_ok.update(dic_var_levs_exclude)

    for k in range(k_inicial, k_final + 1):
        for col in dic_var_levs_exclude_ok.keys():
            if k in dic_var_levs_exclude_ok[col]:
                df_all_cols[col].loc[df_all_cols['k'] == k] = dic_defatult_values[col]

            # df_all_cols_tmp = df_all_cols[['k', col]].copy()
            # if k in dic_var_levs_exclude_ok[col]:
            #     df_all_cols_tmp[col] = pd.Series(dic_defatult_values[col])
            # # else:
            # #     df_all_cols_tmp[col] = pd.Series(df_all_cols[col])
            # df_all_cols = pd.concat([df_all_cols, df_all_cols_tmp], ignore_index=True)

    return df_all_cols


def get_all_minmax_values():

    if k_max == 64:

        all_minmax_values = {'k': [-11.600000000000001, 76.6], 'sl': [-0.19928322752178643, 1.1967689676659643],
         'Tc': [18.113996373799985, 373.6123756357], 'qv': [-0.0035332086576800007, 0.02119925195308],
         'qc': [-0.00019144478858000002, 0.00114866873148], 'qr': [-0.000289795443728, 0.0017387726623679998],
         'qi': [-2.4782610170000003e-05, 0.00014869566102], 'qs': [-0.000438658347798, 0.0026319500867879996],
         'qg': [-9.8722767228e-05, 0.0005923366033679999], 'ni': [-238125.60474200002, 1428753.6284520002],
         'ns': [-559478.0950620001, 3356868.5703720003], 'nr': [-923696.847572, 5542181.0854319995],
         'NG': [-7854.751731060001, 47128.51038636], 'NC': [-15902188.53006, 95413131.18035999],
         'tke': [-59.964000000000006, 359.994], 'omega': [-2.2730674705028, 0.9723025949368]}

        # TODO
         # 'LSRAIN': [-5.45175367804e-05, 0.0003271052206824], 'LSSNOW': [-1.684981073452e-11, 1.0109886440712e-10]}

    elif k_max == 28:

        # TODO FOR TKE
        # 2014 + 2015 - MLP + LSTM  - borde 0
        # {'k': [1, 28], 'sl': [0.00152250216953, 0.994963901883], 'Tc': [63.6765959722, 311.233042513],
        #  'qv': [1.53669990465e-06, 0.0173136957911], 'qc': [0.0, 0.000740654592221], 'qr': [0.0, 0.00080807454392],
        #  'qi': [0.0, 0.000116620887624], 'qs': [0.0, 0.00185267055718], 'qg': [0.0, 0.000135966743237],
        #  'ni': [0.0, 1158020.15791], 'ns': [0.0, 2973056.2935], 'nr': [0.0, 2951432.68921], 'NG': [0.0, 5779.33190043],
        #  'NC': [0.0, 57278488.654], 'kzh': [0.03, 300.0], 'omega': [-1.80612289269, 0.504465053754]}

        # 2014 + 2015 - MLP + LSTM
        all_minmax_values = {'k': [-4.4, 33.4],
                             'sl': [-0.197165777773164, 1.193652181825694],
                             'Tc': [14.165306664039996, 360.74433182116],
                             'qv': [-0.0034608951183344194, 0.020776127609339067],
                             'qc': [-0.0001481309184442, 0.0008887855106652],
                             'qr': [-0.000161614908784, 0.000969689452704],
                             'qi': [-2.33241775248e-05, 0.0001399450651488],
                             'qs': [-0.000370534111436, 0.0022232046686159997],
                             'qg': [-2.71933486474e-05, 0.0001631600918844],
                             'ni': [-231604.03158200003, 1389624.189492],
                             'ns': [-594611.2587, 3567667.5522000003],
                             'nr': [-590286.537842, 3541719.227052],
                             'NG': [-1155.8663800860002, 6935.1982805160005],
                             'NC': [-11455697.730800001, 68734186.3848],
                             'kzh': [-59.964000000000006, 359.994],
                             'omega': [-2.2682404819788, 0.9665826430428],
                             'LSRAIN': [-4.5165380572800005e-05, 0.00027099228343680003],
                             'LSSNOW': [-1.5685887775440002e-10, 9.411532665264e-10]}

    # elif k_max == 18:
    #     all_minmax_values = {'k': [-2.4000000000000004, 21.4], 'si': [-0.16842105263156, 1.19473684210526],
    #                          'Tc': [121.66744111119999, 343.79478323880005],
    #                          'qv': [-0.00357206037104, 0.021459618476240002],
    #                          'qc': [-0.00045647599200200004, 0.002738855952012],
    #                          'qr': [-0.00039292703813800003, 0.002357562228828],
    #                          'qi': [-2.46881644156e-05, 0.0001481289864936],
    #                          'qs': [-0.0004968071975520001, 0.0029808431853120005],
    #                          'qg': [-0.00015900436603420002, 0.0009540261962052],
    #                          'ni': [-288120.643416, 1728723.860496],
    #                          'ns': [-762211.516934, 4573269.101604],
    #                          'nr': [-1225400.800416, 7352404.802496],
    #                          'NG': [-18024.268847020005, 108145.61308212002],
    #                          'NC': [-12525940.75164, 75155644.50984],
    #                          'tke': [-1.2000000000000002, 7.2],
    #                          'kzh': [-59.964000000000006, 359.994],
    #                          'omega': [-2.2562946315226, 0.9304536572056001],
    #                          'EFFCS': [-5.8909090909, 42.3454545454],
    #                          'EFFIS': [-10.400000000000002, 153.4],
    #                          'LSRAIN': [-0.0001374135668268, 0.0008244814009608],
    #                          'LSSNOW': [-2.41726074904e-07, 1.4503564494240001e-06]}

    return all_minmax_values


# def preprocess_lstm2d(X, TimeSteps, stride=1):
#
#     X_samples = list()
#     NumerOfRows = len(X)
#     time_stride = int(TimeSteps * stride)
#     print('time_stride, NumerOfRows', time_stride, NumerOfRows)
#     for i in range(time_stride, NumerOfRows, 1):
#         x_sample = X.iloc[i - time_stride:i:stride]
#         X_samples.append(x_sample)
#
#     X_data = np.array(X_samples)
#     print("X_data.shape=", X_data.shape)
#     return X_data


def get_arr_input_windowed(b, window):
    w = np.empty(((b.shape[0] - window + 1) * window, b.shape[1]))

    for i in range(window, b.shape[0] + 1):
        for ker in range(window):
            w[(i - window) * window + ker, :] = b[i - window + ker, :]
    return w



# Begin program
#

# parameters ============================
modelType = 'MLP_LSTM_LS'
debug_one_line_mode = False
k_max = 64
is_col_k = False
is_PCA = True
timesteps = 5  # use same as bam1d, even with MLP_LSTM_LS or MLP_colk_LSTM_LS with 1 timestep


pd.set_option("display.max_rows", None, "display.max_columns", None)
# evita o erro A value is trying to be set on a copy of a slice from a DataFrame em
pd.options.mode.chained_assignment = None
if is_col_k:
    lev_str = 'Lev'
else:
    lev_str = 'NoLev'

all_minmax_values = get_all_minmax_values()
all_out_columns = ["k", "Tc", "qv", "qc", "qr", "qi", "qs", "qg", "ni", "ns", "nr", "NG", "NC", "LSRAIN", "LSSNOW"]

# timesteps = 1
#
# modelType = 'MLP'

#
# modelType = 'MLP_T'
# modelType = 'CNN_T'
# modelType = 'LSTM'
# lines_per_exec = (k_max * timesteps)
# MODEL MLP_NOLS + LSTM_LS

lines_per_exec = (k_max * timesteps)

s1 = time.time()

model_MLP = None
model_LSTM = None
model = None
all_LSTM_columns = None
all_MLP_columns_inp = None
all_MLP_columns_out = None

if modelType == 'MLP_LSTM_LS' or modelType == 'MLP_colk_LSTM_LS':
    model_MLP = load_model(f'./MLP_LS.h5')
    print(model_MLP.summary())
    model_LSTM = load_model(f'./LSTM_LS.h5')
    print(model_LSTM.summary())
    if is_PCA:
        # later reload the pickle file
        pca_reload = pk.load(open("pca.pkl", 'rb'))

    all_MLP_columns_inp = ["Tc", "qv", "qc", "qr", "qi", "qs", "qg", "ni", "ns", "nr", "NG", "NC", "sl", "tke", "omega"]
    all_MLP_columns_out = ["Tc", "qv", "qc", "qr", "qi", "qs", "qg", "ni", "ns", "nr", "NG", "NC"]
    all_LSTM_columns = ["LSRAIN", "LSSNOW"]
else:  # 'MLP_colk', 'MLP_T', 'CNN_T', 'LSTM'
    model = load_model(f'./{modelType}.h5')
    print(model.summary())

s2 = time.time()
print("Time to load model: ", s2 - s1)

bam1d_path = '/home/denis/_COM_BACKUP/NN_BAM1d/bam1d/work/model/exec'
# bam1d_path = '.'
in_filename = bam1d_path + '/BAM1D_HUMO_in.csv'
out_filename = bam1d_path + '/BAM1D_HUMO_NN_out.csv'

dic_var_levels_exclude_examples, dic_var_levels_exclude_targets = get_dic_levels_excluded(lev_str)

if not debug_one_line_mode:
    try:
        os.remove(in_filename)
    except Exception as e:
        pass
    try:
        os.remove(out_filename)
    except Exception as e:
        pass

# TODO - while not ended BAM1D
while True:
    # s2 = time.time()
    df_input_bam1d = None
    notok = True
    while notok:
        # read all inputs from BAM1D, normalize and remove unecessary cols
        try:
            # time.sleep(1)
            num_lines = sum(1 for line in open(in_filename))
            # time.sleep(0.01)
            # se existe arquivo e tem linhas suficientes
            if num_lines > lines_per_exec:  #TODO deveria ser == lines_per_exec + 1 ?
                df_input_bam1d = pd.read_csv(in_filename, sep=",")
                print('shape=', df_input_bam1d.shape)
                notok = df_input_bam1d.shape[0] < lines_per_exec

            # TODO - se for usar para comparar, ler do output total
            # df_output_bam1d = pd.read_csv('./BAM1D_HUMO_out.csv', sep=","
        except Exception as e:
            pass
            # print('waiting for BAM1D write input and output files')

    if not debug_one_line_mode:
        os.remove(in_filename)
    df_input_bam1d, all_input_cols = preprocess_features(df_input_bam1d)
    # df_output_bam1d = preprocess_targets(df_output_bam1d)

    # print('Input FROM BAM \n', df_input_bam1d.shape)
    print('Input denorm \n', df_input_bam1d)
    df_input_bam1d_norm = normalize_linear_scale(df_input_bam1d, all_minmax_values)
    # print('Input norm \n', df_input_bam1d_norm)

    if modelType == 'MLP_colk':
        if is_col_k:
            df_input_bam1d_norm = get_df_col_k(df_input_bam1d_norm, 1, k_max, dic_var_levels_exclude_examples)
            # print('Input norm exclude levs \n', df_input_bam1d_norm.shape)
        np_input_bam1d_norm = df_input_bam1d_norm.to_numpy()

    elif modelType == 'MLP_colk_LSTM_LS' or modelType == 'MLP_LSTM_LS':

        if is_col_k:
            df_input_bam1d_norm_MLP = df_input_bam1d_norm[['k'] + all_MLP_columns_inp]
            df_input_bam1d_norm_MLP = get_df_col_k(df_input_bam1d_norm_MLP, 1, k_max, dic_var_levels_exclude_examples)
        else:
            df_input_bam1d_norm_MLP = df_input_bam1d_norm[all_MLP_columns_inp]

        if is_PCA:
            np_input_bam1d_norm_MLP = pca_reload.transform(df_input_bam1d_norm_MLP)
        else:
            np_input_bam1d_norm_MLP = df_input_bam1d_norm_MLP.to_numpy()

        df_input_bam1d_norm_LSTM_temp = df_input_bam1d_norm[all_LSTM_columns]
        df_input_bam1d_norm_LSTM = pd.DataFrame(columns=all_LSTM_columns)
        # value repeated ecah k
        for w in range(timesteps):
            df_input_bam1d_norm_LSTM = df_input_bam1d_norm_LSTM.append(df_input_bam1d_norm_LSTM_temp.iloc[(w * k_max) + k_max - 1, :], ignore_index=True)

        np_input_bam1d_norm_LSTM_temp = df_input_bam1d_norm_LSTM.to_numpy()
        np_input_bam1d_norm_LSTM = get_arr_input_windowed(np_input_bam1d_norm_LSTM_temp, timesteps)
        new_samples = int(np_input_bam1d_norm_LSTM.shape[0] / timesteps)
        np_input_bam1d_norm_LSTM = np_input_bam1d_norm_LSTM[0:new_samples * timesteps, :]
        np_input_bam1d_norm_LSTM = np_input_bam1d_norm_LSTM.reshape(new_samples, timesteps, np_input_bam1d_norm_LSTM.shape[1])


    else:  # in ['MLP_T', 'CNN_T', 'LSTM'']:
        np_temp = df_input_bam1d_norm.to_numpy()
        np_input_bam1d_norm = get_arr_input_windowed(np_temp, timesteps)
        new_samples = int(np_input_bam1d_norm.shape[0] / timesteps)
        np_input_bam1d_norm = np_input_bam1d_norm[0:new_samples * timesteps, :]
        np_input_bam1d_norm = np_input_bam1d_norm.reshape(new_samples, timesteps, np_input_bam1d_norm.shape[1])

    # Preditction
    #
    if modelType == 'MLP_colk_LSTM_LS' or modelType == 'MLP_LSTM_LS':

        # MLP part
        np_predict_MLP = model_MLP.predict(np_input_bam1d_norm_MLP, use_multiprocessing=True)
        if is_col_k:
            df_predict_bam1d_norm_MLP = get_dffull_from_npcol_k(np_predict_MLP[0], all_MLP_columns_out, 1, k_max,
                                                            dic_var_levs_exclude=dic_var_levels_exclude_targets,
                                                            dic_defatult_values=get_dic_default_values())
        else:
            df_predict_bam1d_norm_MLP = pd.DataFrame(np_predict_MLP, columns=all_MLP_columns_out)

        # print('MLP input norm \n', np_input_bam1d_norm_MLP)
        # print('MLP predict norm \n', np_predict_MLP)

        df_predict_bam1d_MLP = denormalize_linear_scale(df_predict_bam1d_norm_MLP, all_minmax_values, is_col_k=is_col_k)
        # all predicted denormalized has minimun of zero
        df_predict_bam1d_MLP.clip(lower=0.0, inplace=True)

        # LSTM part
        np_predict_LSTM = model_LSTM.predict(np_input_bam1d_norm_LSTM, use_multiprocessing=True)
        print('LSTM input norm \n', np_input_bam1d_norm_LSTM)
        print('LSTM predict norm \n', np_predict_LSTM)
        np_predict_LSTM = np.repeat(np_predict_LSTM, k_max, axis=0)
        df_predict_bam1d_norm_LSTM = pd.DataFrame(np_predict_LSTM, columns=all_LSTM_columns)
        df_predict_bam1d_LSTM = denormalize_linear_scale(df_predict_bam1d_norm_LSTM, all_minmax_values,
                                                              is_col_k=is_col_k)

        df_predict_bam1d = df_predict_bam1d_MLP.join(df_predict_bam1d_LSTM)
        print('MLP + LS prediction denorm \n', df_predict_bam1d)
        df_predict_bam1d.to_csv(out_filename, header=False, index=False)

        # print('Dataframe actual \n', df_output_bam1d)
        # print('MAE \n', mae(df_output_bam1d, df_predict_bam1d_norm))
        # print('BIAS \n', bias(df_output_bam1d, df_predict_bam1d_norm))
        # print('RMSE \n', rmse(df_output_bam1d, df_predict_bam1d_norm))

    else:
        if is_col_k == True:
            np_col_k_predict = model.predict(np_input_bam1d_norm, use_multiprocessing=True)
            df_predict_bam1d_norm = get_dffull_from_npcol_k(np_col_k_predict[0], all_out_columns, 1, k_max,
                                                            dic_var_levs_exclude=dic_var_levels_exclude_targets,
                                                            dic_defatult_values=get_dic_default_values())

        else:  # TODO
            # df_predict_bam1d_norm = ?
            pass

        df_predict_bam1d_norm = denormalize_linear_scale(df_predict_bam1d_norm, all_minmax_values, is_col_k=is_col_k)

        # print('Predicted \n', np_col_k_predict)
        # s4 = time.time()
        # print("Time to Prediction: ", s4-s3)
        # create df default values and read values from prediction
        # print('Post default values \n', df_predict_bam1d_norm)
        # s5 = time.time()
        # print("Time to create denorm df: ", s5 - s4)
        print('Dataframe prediction denorm \n', df_predict_bam1d_norm)

        df_predict_bam1d_norm = set_default_values(df_predict_bam1d_norm, 1, k_max,
                                                   dic_var_levs_exclude=dic_var_levels_exclude_targets,
                                                   dic_defatult_values=get_dic_default_values())
        print('Dataframe prediction denorm default values \n', df_predict_bam1d_norm)

        df_predict_bam1d_norm.to_csv(out_filename, header=False, index=False)

        # print('Dataframe actual \n', df_output_bam1d)
        # print('MAE \n', mae(df_output_bam1d, df_predict_bam1d_norm))
        # print('BIAS \n', bias(df_output_bam1d, df_predict_bam1d_norm))
        # print('RMSE \n', rmse(df_output_bam1d, df_predict_bam1d_norm))
