import pandas as pd
import scipy as sp
import numpy as np
import xarray as xa
from matplotlib import pyplot as plt

import pyro
import pyro.distributions as dist
import pyro.distributions.constraints as constraints
from pyro.optim import Adam
from pyro.infer import SVI, TraceEnum_ELBO, config_enumerate, infer_discrete
import torch

from tqdm import tqdm
import scipy as sp
from wisdom.data import load_responses

def prepare_data(data):
    woc_stacked = data.stack().reset_index()
    cat_enm = woc_stacked["ENM"].astype("category")
    cat_hazard = woc_stacked["hazard"].astype("category")
    cat_expert = woc_stacked["expert"].astype("category")
    y = woc_stacked[0].astype(np.float32)

    y_t = torch.from_numpy(y.values)
    cat_enm_t = torch.from_numpy(cat_enm.cat.codes.values.copy()).long()
    cat_hazard_t = torch.from_numpy(cat_hazard.cat.codes.values.copy()).long()
    cat_expert_t = torch.from_numpy(cat_expert.cat.codes.values.copy()).long()

    woc_shape = (cat_expert_t.max().item() + 1, cat_hazard_t.max().item() + 1, cat_enm_t.max().item()  + 1)
    woc_t = torch.full(woc_shape, float('nan'))
    woc_t[cat_expert_t, cat_hazard_t, cat_enm_t] = y_t
    woc_mask_t = ~woc_t.isnan()
    woc_filled_t = woc_t.nan_to_num(0)
    return woc_filled_t, woc_mask_t, cat_enm.cat.categories, cat_hazard.cat.categories, cat_expert.cat.categories

@config_enumerate
def model(y=None, mask=None):
    n_experts, n_hazards, n_enms = y.shape
    enm_plate = pyro.plate("enm", size=n_enms, dim=-1)
    hazard_plate = pyro.plate("hazard", size=n_hazards, dim=-2)
    expert_plate = pyro.plate("expert", size=n_experts, dim=-3)
    with expert_plate:
        expertise = pyro.sample("expertise", dist.Normal(1, 1))
    
    with hazard_plate:
        difficulty = torch.exp(pyro.sample("difficulty", dist.Normal(1, 1)))
    
    with enm_plate:
        p = pyro.sample("p", dist.Beta(1, 1))
        with hazard_plate:
            hidden_label = pyro.sample("hidden_label", dist.Bernoulli(probs=p))
            with expert_plate:
                sigma = torch.sigmoid(expertise * difficulty)
                obs_probs = (sigma**hidden_label)*((1-sigma)**(1-hidden_label))
                with pyro.poutine.mask(mask=mask):
                    pyro.sample("obs", dist.Bernoulli(probs=obs_probs), obs=y)

def guide(y=None, mask=None):
    n_experts, n_hazards, n_enms = y.shape
    p_a = pyro.param("p_a", torch.ones(size=(n_enms, )), constraint=constraints.positive)
    p_b = pyro.param("p_b", torch.ones(size=(n_enms, )), constraint=constraints.positive)
    expert_loc = pyro.param("expert_loc", torch.ones(size=(n_experts, 1, 1)))
    expert_scale = pyro.param("expert_scale", torch.ones(size=(n_experts, 1, 1)), constraint=constraints.positive)
    hazard_loc = pyro.param("hazard_loc", torch.ones(size=(n_hazards, 1)))
    hazard_scale = pyro.param("hazard_scale", torch.ones(size=(n_hazards, 1)), constraint=constraints.positive)
    enm_plate = pyro.plate("enm", size=n_enms, dim=-1)
    hazard_plate = pyro.plate("hazard", size=n_hazards, dim=-2)
    expert_plate = pyro.plate("expert", size=n_experts, dim=-3)
    with expert_plate:
        expertise = pyro.sample("expertise", dist.Normal(expert_loc, expert_scale))
    with hazard_plate:
        difficulty = pyro.sample("difficulty", dist.Normal(hazard_loc, hazard_scale))
    with enm_plate:
        p = pyro.sample("p", dist.Beta(p_a, p_b))

def fit_responses(y_t, mask_t, n_steps=10_000):
    pyro.clear_param_store()
    pyro.enable_validation(True)

    # setup the optimizer
    adam_params = {"lr": 1e-3, "betas": (0.9, 0.999)}
    optimizer = Adam(adam_params)
    # setup the inference algorithm
    svi = SVI(model, guide, optimizer, loss=TraceEnum_ELBO(max_plate_nesting=3))
    # Initializes param store.
    svi.loss( model, guide, y_t, mask_t)  

    losses = []
    pbar = tqdm(range(n_steps))
    # fit the model
    for i in pbar:
        loss = svi.step(y_t, mask_t)
        losses.append(loss)
        if i % 100 == 0:
            pbar.set_postfix_str("ELBO: {:.2f}".format(loss))
    return losses

def sample_response_model(n_samples=1000):
    # sample hidden concern labels
    sample_predict = infer_discrete(model, first_available_dim=-4)
    preds = []
    for _ in range(n_samples):
        preds.append(pyro.poutine.trace(sample_predict).get_trace(y_t, mask_t).nodes["hidden_label"]["value"].numpy())
    mode_z, uncertainty_z = sp.stats.mode(np.array(preds).reshape((1000, -1)), axis=0)
    
    # sample posterior distributions
    p_s = []
    expertise_s = []
    difficulty_s = []
    for i in range(n_samples):
        sample = pyro.poutine.trace(guide).get_trace(y_t, mask_t)
        p_s.append(sample.nodes["p"]["value"].detach().numpy())
        expertise_s.append(sample.nodes["expertise"]["value"].detach().numpy())
        difficulty_s.append(sample.nodes["difficulty"]["value"].detach().numpy())
    
    samples = pd.DataFrame(p_s, columns=cat_enm).T.stack().reset_index()
    expert_scores = pd.DataFrame(np.array(expertise_s).squeeze(), columns=cat_expert).T.stack().reset_index()
    difficulty_scores = pd.DataFrame(np.exp(np.array(difficulty_s).squeeze()), columns=cat_hazard).T.stack().reset_index()
    return preds, mode_z, uncertainty_z, samples, expert_scores, difficulty_scores

if __name__ == "__main__":
    woc = load_responses("./data/")
    woc = woc[woc.index.get_level_values(1) != "Overall.hazardous"]

    enms = woc.index.get_level_values(0).unique().to_list()
    hazards = woc.index.get_level_values(1).unique().to_list()
    experts = woc.columns.to_list()

    y_t, mask_t, cat_enm, cat_hazard, cat_expert = prepare_data(woc)

    losses = fit_responses(y_t, mask_t)
    fig, ax = plt.subplots(figsize=(16, 9))
    ax.plot(losses)
    ax.set_xlabel("step")
    ax.set_ylabel("loss")
    fig.savefig("outputs/response_model_fit.png")

    preds, mode_z, uncertainty_z, samples, expert_scores, difficulty_scores = sample_response_model()

    
    # get posterior parameters
    concerns_a = pyro.get_param_store().get_param("p_a").detach().numpy()
    concerns_b = pyro.get_param_store().get_param("p_b").detach().numpy()
    expert_loc = pyro.get_param_store().get_param("expert_loc").detach().numpy().squeeze()
    expert_scale = pyro.get_param_store().get_param("expert_scale").detach().numpy().squeeze()
    hazard_loc = pyro.get_param_store().get_param("hazard_loc").detach().numpy().squeeze()
    hazard_scale = pyro.get_param_store().get_param("hazard_scale").detach().numpy().squeeze()
    
    # store results
    hidden_labels = xa.Dataset(
        {
            "samples": (("draw", "hazard", "enm"), np.array(preds)),
            "inferred_labels": (("hazard", "enm"), mode_z.reshape((17, 134))),
            "label_uncertainty": (("hazard", "enm"), uncertainty_z.reshape((17, 134))),
            "posterior_global_concerns": (("draw", "enm"), samples.pivot(columns="level_0", index="level_1", values=0)),
            "posterior_expert_scores": (("draw", "expert"), expert_scores.pivot(columns="level_0", index="level_1", values=0)),
            "posterior_difficulty_score": (("draw", "hazard"), difficulty_scores.pivot(columns="level_0", index="level_1", values=0)),
            "posterior_concerns_a": (("enm",), concerns_a),
            "posterior_concerns_b": (("enm",), concerns_b),
            "posterior_expert_loc": (("expert",), expert_loc),
            "posterior_expert_scale": (("expert",), expert_scale),
            "posterior_hazard_loc": (("hazard",), hazard_loc),
            "posterior_hazard_scale": (("hazard",), hazard_scale),
        },
        coords = {"hazard": cat_hazard, "enm": cat_enm, "expert": cat_expert}
    )

    hidden_labels.to_netcdf("./outputs/model_inference_anon.nc")