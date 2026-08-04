# Module for a lightweight RNN model class only for prediction

import copy
import warnings
import numpy as np
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
        inputs, outputs, state_specs = self._build_model(params)
        super().__init__(inputs=inputs, outputs=outputs, **kwargs)

        self.params = params
        self.state_specs = state_specs
        self._cycle_states = None

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
        initial_state_inputs = []
        final_state_outputs = []
        state_specs = []
        for i, layer_type in enumerate(params["hidden_layers"]):
            units = params["hidden_units"][i]
            activation = params["hidden_activation"][i]

            if layer_type == "dense":
                x = layers.Dense(units=units, activation=activation)(x)
            elif layer_type == "dropout":
                x = layers.Dropout(params["dropout"])(x)
            elif layer_type == "rnn":
                h0 = Input(shape=(units,), name=f"rnn_{i}_h0")
                x, h = layers.SimpleRNN(
                    units=units,
                    activation=activation,
                    dropout=params["dropout"],
                    recurrent_dropout=params["recurrent_dropout"],
                    stateful=False,
                    return_sequences=True,
                    return_state=True
                )(x, initial_state=[h0])
                initial_state_inputs.append(h0)
                final_state_outputs.append(h)
                state_specs.append({
                    "layer_index": i,
                    "layer_type": "rnn",
                    "state_names": ["h"],
                    "units": units,
                })                
            elif layer_type == "lstm":
                h0 = Input(shape=(units,), name=f"lstm_{i}_h0")
                c0 = Input(shape=(units,), name=f"lstm_{i}_c0")
                x, h, c = layers.LSTM(
                    units=units,
                    activation=activation,
                    dropout=params["dropout"],
                    recurrent_dropout=params["recurrent_dropout"],
                    stateful=False,
                    return_sequences=True,
                    return_state=True
                )(x, initial_state=[h0, c0])
                initial_state_inputs.extend([h0, c0])
                final_state_outputs.extend([h, c])
                state_specs.append({
                    "layer_index": i,
                    "layer_type": "lstm",
                    "state_names": ["h", "c"],
                    "units": units,
                })                
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

        return x, initial_state_inputs, final_state_outputs, state_specs

    @classmethod
    def _build_model(cls, params):
        """
        Build a flexible sequence-to-sequence prediction graph.
        """
        inputs = Input(batch_shape=(None, None, params["n_features"]))
        x, initial_state_inputs, final_state_outputs, state_specs = cls._build_hidden_layers(inputs, params)

        if params["output_layer"] == "dense":
            predictions = layers.Dense(
                units=params["output_dimension"],
                activation=params["output_activation"],
            )(x)
        else:
            raise ValueError(f"Unsupported output layer type: {params['output_layer']}")

        return [inputs] + initial_state_inputs, [predictions] + final_state_outputs, state_specs

    @classmethod
    def from_weights(cls, params, weights_path):
        model = cls(params=params)
        model.load_weights(weights_path)
        return model

    def _zero_cycle_states(self, batch_size, dtype=np.float32):
        """
        Create zero recurrent states matching the model's recurrent layers.
        """
        states = []
        for spec in self.state_specs:
            for _ in spec["state_names"]:
                states.append(np.zeros((batch_size, spec["units"]), dtype=dtype))
        return states

    def _validate_cycle_states(self, states):
        """
        Validate and flatten recurrent states supplied to predict_cycle.
        """
        if states is None:
            return None

        states = list(states)
        n_expected = sum(len(spec["state_names"]) for spec in self.state_specs)
        if len(states) != n_expected:
            raise ValueError(f"Expected {n_expected} recurrent state arrays, got {len(states)}.")
        return states

    def reset_cycle_states(self):
        """
        Clear stored recurrent states.
        """
        self._cycle_states = None

    def predict_cycle(self, X, reset_state=False, initial_states=None, return_states=False, **kwargs):
        """
        Stores recurrent states after prediction and continues from stored states if they exist. Used for operational prediction where input data might come in cycles

        Args
        =========
        X: ndarray, input data (nbatch, ntime, nfeatures)
        reset_state: bool, whether to reset recurrent states (to zeros by default). Use if predicting at a new location or time
        initial_states: list, optional flat list of recurrent states. If None,
            use stored states when available.
        return_states: bool, whether to return final recurrent states as a flat
            list ordered by recurrent layer, with [h] for SimpleRNN and [h, c]
            for LSTM.
        """
        if initial_states is not None:
            cycle_states = self._validate_cycle_states(initial_states)
        elif reset_state or self._cycle_states is None:
            x_array = np.asarray(X)
            cycle_states = self._zero_cycle_states(batch_size=x_array.shape[0], dtype=x_array.dtype)
        else:
            cycle_states = self._validate_cycle_states(self._cycle_states)

        outputs = super().predict([X] + cycle_states, **kwargs)
        if not isinstance(outputs, (list, tuple)):
            outputs = [outputs]

        predictions = outputs[0]
        self._cycle_states = outputs[1:]

        if return_states:
            return predictions, self._cycle_states
        return predictions
        

def warp_weights(weights0, bi_warp, bf_warp):
    """
    Given LSTM layer weights and time-warp parameters, return a new list
    of time-warped LSTM weights without modifying the input weights.
    """
    # Copy all arrays to avoid mutating the originals
    w_warped = [w.copy() for w in weights0]
    # Bias vector (Keras LSTM layout: [i, f, c, o])
    b = w_warped[2]
    # Infer number of LSTM units from bias length
    if b.ndim != 1 or b.shape[0] % 4 != 0:
        raise ValueError("Unexpected LSTM bias shape.")
    lstm_units = b.shape[0] // 4
    # Input gate biases (i)
    b[0:lstm_units] += bi_warp
    # Forget gate biases (f)
    b[lstm_units:2 * lstm_units] += bf_warp

    return w_warped
