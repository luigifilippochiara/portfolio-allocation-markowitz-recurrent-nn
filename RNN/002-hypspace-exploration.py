import h5py
from keras import layers, optimizers
from keras.layers import Input, Dense, LSTM, Dropout,Activation
from keras.models import Model, Sequential
from keras.callbacks import EarlyStopping, ModelCheckpoint
import numpy as np
import telegram_send
from hyperopt import fmin, tpe, hp, STATUS_OK, Trials
import json
from keras import backend as K

baseDataset = 'month-month-month-ret'

# check nr of gpus available
gpus = K.tensorflow_backend._get_available_gpus()
print('Available gpus', gpus)

# import data
def returnData(windowSize):
    with h5py.File('datasets-{}/train-{}.h5'.format(baseDataset, windowSize), 'r') as hf:
        xTrain = hf['x'][:]
        yTrain = hf['y'][:]
    with h5py.File('datasets-{}/val-{}.h5'.format(baseDataset, windowSize), 'r') as hf:
        xVal = hf['x'][:]
        yVal = hf['y'][:]
    with h5py.File('datasets-{}/test-{}.h5'.format(baseDataset, windowSize), 'r') as hf:
        xTest = hf['x'][:]
        yTest = hf['y'][:]
    return xTrain, yTrain, xVal, yVal, xTest, yTest

# hyperparameters space
space = {
    'neurons': hp.choice('neurons', [32,64,128,256]),
    'dropout': hp.choice('dropout', [0,0.1,0.2,0.3,0.4]),
    'batchSize': hp.choice('batchSize', [32,64]),
    'twoLayers': hp.choice('twoLayers', [True,False]),
    'windowSize': hp.choice('windowSize', np.arange(9,61))
}

currentEval = 0
maxEvaluations = 1000
# telegram notification steps
notificationSteps = np.round(np.quantile(np.arange(maxEvaluations), np.linspace(0,1,20)))

# model definition
def create_model(space):
    global currentEval
    global maxEvaluations

    # clear memory to avoid "out of memory" error
    if K.backend() == 'tensorflow':
        K.clear_session()
    
    # load data
    xTrain, yTrain, xVal, yVal, xTest, yTest = returnData(space['windowSize'])
    
    # define model
    model = Sequential()
    if space['twoLayers']:
        model.add(LSTM(space['neurons'], return_sequences=True, input_shape=(xTrain.shape[1], xTrain.shape[2]), activation='tanh'))
        model.add(Dropout(space['dropout']))
        model.add(LSTM(space['neurons'], activation='tanh'))
    else:
        model.add(LSTM(space['neurons'], return_sequences=False, input_shape=(xTrain.shape[1], xTrain.shape[2]), activation='tanh'))
        model.add(Dropout(space['dropout']))
    model.add(Dense(units=xTrain.shape[2]))
    model.compile(loss='mse', optimizer='adam', metrics=['mae'])

    
    # callbacks
    # stop if loss does not decrease after x epochs
    callbacks = [EarlyStopping(patience=10, monitor='val_loss', min_delta=0, mode='min')]

    result = model.fit(xTrain, yTrain,
              batch_size=space['batchSize'],
              epochs=200,
              verbose=2,
              validation_data=(xVal, yVal),
              callbacks=callbacks,
              shuffle=False)
    
    #get the highest validation accuracy of the training epochs
    validation_loss = np.min(result.history['val_loss']) 
    print('Best validation loss of evaluation #{}/{}: {}'.format(currentEval, maxEvaluations, validation_loss))
    print()
    
    # notify via telegram every 25% of the process
    if currentEval in notificationSteps:
        telegram_send.send(['Hyperparameter space exploration reached {} of {} evaluations'.format(currentEval, maxEvaluations)])
    
    currentEval += 1
    return {'loss': validation_loss, 'status': STATUS_OK}

# start hyperparameter space exploration
trials = Trials()
best = fmin(create_model, space, algo=tpe.suggest, max_evals=maxEvaluations, trials=trials)

# store results
results = {
    'best': best,
    'trials': trials.trials,
    'results': trials.results,
    'best_trial': trials.best_trial
}
with open('results-{}.json'.format(baseDataset), 'w') as fp:
    json.dump(results, fp, default=str)

# send finish notification
telegram_send.send(['Hyperparameter search finished'])