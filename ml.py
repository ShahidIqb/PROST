import joblib

def ens_predict(Xtest):
    """Load stored models and predict the outputs from input features"""
    extraTrees_model = joblib.load('models/extraTrees_forward_reverse')
    xgboost_model = joblib.load('models/xgb_forward_reverse')
    #lin_model = joblib.load('linear_forward_reverse')
    extr_pred_y = extraTrees_model.predict(Xtest)
    xgb_pred_y = xgboost_model.predict(Xtest)
    #comb_extr_xgb = numpy.c_[extr_pred_y, xgb_pred_y]
    #ens_pred_y = lin_model.predict(comb_extr_xgb)
    ens_wtAvg = 0.15*extr_pred_y + 0.85 * xgb_pred_y
    return ens_wtAvg