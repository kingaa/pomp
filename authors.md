---
title: pomp authors
layout: pomp
---

## Authors

{% for author in site.data.authors %}> [{{ author.name }}]({{ author.url }})  
{% endfor %}
