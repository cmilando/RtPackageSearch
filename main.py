import requests


def by_search_term(query):
    per_page = 100
    page = 1
    url = f"https://api.github.com/search/repositories?q={query}&per_page={per_page}&page={page}"

    response = requests.get(url)
    data = response.json()

    for repo in data.get("items", []):
        print(f"{query}| {repo['full_name']}| {repo['html_url']} | {repo['description']} | {repo['updated_at']}")

    page = 2
    url = f"https://api.github.com/search/repositories?q={query}&per_page={per_page}&page={page}"

    response = requests.get(url)
    data = response.json()

    for repo in data.get("items", []):
        print(f"{query}| {repo['full_name']}| {repo['html_url']} | {repo['description']} | {repo['updated_at']}")


if __name__ == "__main__":

    by_search_term("reproduction numbers")
    by_search_term("reproduction number")

    by_search_term("reproductive numbers")
    by_search_term("reproductive number")

    by_search_term("reproduction rate")
    by_search_term("reproduction rates")



