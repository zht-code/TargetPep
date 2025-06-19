import torch
from torch.utils.data._utils.collate import default_collate

class PaddingCollate:
    def __init__(self, pad_keys, pad_value=0.0):
        """
        Initialize the PaddingCollate instance with specified keys to pad and a padding value.
        
        Args:
        - pad_keys (list of str): List of keys that need padding.
        - pad_value (float): The value used for padding numeric data. Default is 0.0.
        """
        self.pad_keys = pad_keys
        self.pad_value = pad_value

    def pad_tensor(self, tensor, max_len):
        """
        Pad the tensor to the maximum length in the batch.

        Args:
        - tensor (torch.Tensor): The tensor to pad.
        - max_len (int): The maximum length to pad to.

        Returns:
        - torch.Tensor: The padded tensor.
        """
        # if tensor.size(0) < max_len:
        #     padding_size = (max_len - tensor.size(0),) + tuple(tensor.shape[1:])
        #     padding = torch.full(padding_size, self.pad_value, dtype=tensor.dtype, device=tensor.device)
        #     tensor = torch.cat([tensor, padding], dim=0)
        if tensor.dim() == 1:  # Simple 1D padding
            if tensor.size(0) < max_len:
                padding_size = (max_len - tensor.size(0),)
                padding = torch.full(padding_size, self.pad_value, dtype=tensor.dtype, device=tensor.device)
                tensor = torch.cat([tensor, padding], dim=0)
        else:  # Multidimensional padding
            padding_size = [0] * (2 * len(tensor.shape))
            padding_size[-1] = max_len - tensor.size(0)  # Pad the first dimension
            tensor = torch.nn.functional.pad(tensor, padding_size, value=self.pad_value)
        return tensor

    def collate(self, batch):
        """
        Collate the batch, padding specified tensors to the same size.

        Args:
        - batch (list): A list of dictionaries containing the data samples.

        Returns:
        - dict: A dictionary containing the collated data with all specified tensors padded to the same size.
        """
        # Find the maximum length for each tensor key that needs padding
        max_lens = {key: max(x[key].size(0) for x in batch if key in x) for key in self.pad_keys}

        # Pad each tensor that needs padding
        for data in batch:
            for key in self.pad_keys:
                if key in data and data[key].size(0) < max_lens[key]:
                    data[key] = self.pad_tensor(data[key], max_lens[key])
        # # Ensure all tensors for each key have the same shape
        # uniform_shapes = {}
        # for key in self.pad_keys:
        #     shapes = {data[key].shape for data in batch if key in data}
        #     if len(shapes) > 1:
        #         max_shape = tuple(max(sizes) for sizes in zip(*shapes))
        #         for data in batch:
        #             if key in data:
        #                 data[key] = self.pad_tensor(data[key], max_shape[0])           

        return default_collate(batch)

    def __call__(self, batch):
        return self.collate(batch)
