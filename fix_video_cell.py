import json

notebook_path = 'notebooks/party_cartograms_zweitstimme_prozent_nonglobal_density_scaling.ipynb'

with open(notebook_path) as f:
    nb = json.load(f)

# Remove the last cell (the broken video cell we added)
if nb['cells'] and '## 8. Result' in ''.join(nb['cells'][-1].get('source', [])):
    nb['cells'].pop()

# Add a raw HTML cell for the video (raw cells are passed through by nbsphinx)
video_cell = {
    "cell_type": "raw",
    "metadata": {
        "raw_mimetype": "text/restructuredtext"
    },
    "source": [
        ".. raw:: html\n",
        "\n",
        "   <h2>8. Result: Animated Video</h2>\n",
        "   <p>The animation shows each party's cartogram morphing from the original geography.</p>\n",
        "   <video controls loop autoplay muted width=\"100%\" style=\"max-width: 800px;\">\n",
        "     <source src=\"animation/zweitstimme_parties_nonglobal_scaling.mp4\" type=\"video/mp4\">\n",
        "     Your browser does not support the video tag.\n",
        "   </video>\n"
    ]
}

nb['cells'].append(video_cell)

with open(notebook_path, 'w') as f:
    json.dump(nb, f, indent=1)

print("Fixed video cell in notebook")
