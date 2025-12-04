import json

notebook_path = '../notebooks/03_party_cartograms_zweitstimme_prozent_nonglobal_density_scaling.ipynb'

with open(notebook_path) as f:
    nb = json.load(f)

# Add a markdown cell to display the video
video_cell = {
    "cell_type": "markdown",
    "metadata": {},
    "source": [
        "## 8. Result: Animated Video\n",
        "\n",
        "The animation shows each party's cartogram morphing from the original geography.\n",
        "\n",
        "```{raw} html\n",
        "<video controls width=\"100%\">\n",
        "  <source src=\"animation/zweitstimme_parties_nonglobal_scaling.mp4\" type=\"video/mp4\">\n",
        "  Your browser does not support the video tag.\n",
        "</video>\n",
        "```"
    ]
}

nb['cells'].append(video_cell)

with open(notebook_path, 'w') as f:
    json.dump(nb, f, indent=1)

print("Added video cell to notebook")
