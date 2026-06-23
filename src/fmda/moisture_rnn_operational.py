# Module for a lightweight RNN model class only for prediction

import copy
import warnings
from tensorflow.keras import layers, Model
from tensorflow.keras.layers import Input
from tensorflow.keras.optimizers import Adam


class OperationalRNNPredictor(Model):
    """
    Lightweight RNN model for operational prediction.

    Builds the same flexible architecture as RNN_Flexible from a params dict,
    but omits training callbacks, history plotting, and evaluation helpers.
    """

    def __init__(self, params: dict, compile_model: bool = False, **kwargs):
        params = self._check_params(params)
        inputs, outputs = self._build_model(params)
        super().__init__(inputs=inputs, outputs=outputs, **kwargs)

        self.params = params

        if compile_model:
            optimizer = Adam(learning_rate=self.params["learning_rate"])
            self.compile(loss="mean_squared_error", optimizer=optimizer)

    @staticmethod
    def _check_params(params):
        """
        Force flexible, stateless sequence prediction settings.
        """
        if params is None:
            raise ValueError("params must be provided for OperationalRNNPredictor.")

        params = copy.deepcopy(dict(params))
        params["n_features"] = len(params["features_list"])

        if params.get("timesteps") is not None:
            warnings.warn("timesteps should be None for flexible prediction. Overriding to None.")
            params["timesteps"] = None

        if params.get("return_sequences") is not True:
            warnings.warn("return_sequences should be True for operational prediction. Overriding to True.")
            params["return_sequences"] = True

        if params.get("stateful") is not False:
            warnings.warn("stateful should be False for operational prediction. Overriding to False.")
            params["stateful"] = False

        layer_count = len(params["hidden_layers"])
        if len(params["hidden_units"]) != layer_count or len(params["hidden_activation"]) != layer_count:
            raise ValueError("hidden_layers, hidden_units, and hidden_activation must have the same length.")

        return params

    @staticmethod
    def _build_hidden_layers(x, params):
        """
        Build hidden layers from the parallel params lists.
        """
        for i, layer_type in enumerate(params["hidden_layers"]):
            units = params["hidden_units"][i]
            activation = params["hidden_activation"][i]

            if layer_type == "dense":
                x = layers.Dense(units=units, activation=activation)(x)
            elif layer_type == "dropout":
                x = layers.Dropout(params["dropout"])(x)
            elif layer_type == "rnn":
                x = layers.SimpleRNN(
                    units=units,
                    activation=activation,
                    dropout=params["dropout"],
                    recurrent_dropout=params["recurrent_dropout"],
                    stateful=False,
                    return_sequences=True,
                )(x)
            elif layer_type == "lstm":
                x = layers.LSTM(
                    units=units,
                    activation=activation,
                    dropout=params["dropout"],
                    recurrent_dropout=params["recurrent_dropout"],
                    stateful=False,
                    return_sequences=True,
                )(x)
            elif layer_type == "attention":
                x = layers.Attention()([x, x])
            elif layer_type == "conv1d":
                kernel_size = params.get("kernel_size", 3)
                x = layers.Conv1D(
                    filters=units,
                    kernel_size=kernel_size,
                    activation=activation,
                    padding="same",
                )(x)
            else:
                raise ValueError(f"Unrecognized layer type: {layer_type}")

        return x

    @classmethod
    def _build_model(cls, params):
        """
        Build a flexible sequence-to-sequence prediction graph.
        """
        inputs = Input(batch_shape=(None, None, params["n_features"]))
        x = cls._build_hidden_layers(inputs, params)

        if params["output_layer"] == "dense":
            outputs = layers.Dense(
                units=params["output_dimension"],
                activation=params["output_activation"],
            )(x)
        else:
            raise ValueError(f"Unsupported output layer type: {params['output_layer']}")

        return inputs, outputs

    @classmethod
    def from_weights(cls, params, weights_path):
        model = cls(params=params)
        model.load_weights(weights_path)
        return model


